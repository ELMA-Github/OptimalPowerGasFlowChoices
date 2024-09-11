include("postprocessing.jl")


function SCP(
    ES::EnergySystem,
    config_dict_solver::Dict{Symbol,<:Any}=Dict{Symbol,Any}();
    order::Int64=1,
    k_max::Int64=100, 
    augmentation_method::String, # augmentation norm
    δ_init::Real=10^-2, # Penalty factor augmentation
    δ_update::Real=2, # Update factor for penalty
    δ_max::Real=10^5, # Maximum penalty parameter
    ϕ_max::Real=10^-8, # Maximum violation
    PDE_type::Int=3 # PDE model type
    )
    """
    Sequential convex programming
    Attribute
    k::Iteration number
    # Define SLP parameters
    δ = [0, 10^-4] # Penalty factor augmentation
    δ_update = 2 # Update factor for penalty
    δ_max = 10^5 # Maximum penalty parameter
    PDE_type = {1,2,3} # PDE model type
    """

    M_base, Π_base, Cost_base, Γ_base =
        ES.bases_dict[:M_base], ES.bases_dict[:Π_base], ES.bases_dict[:Cost_base], ES.bases_dict[:Γ_base]
    # δ_init, δ_max = M_base/Cost_base * [δ_init, δ_max]

    ϕ_max = ϕ_max/Γ_base
    @info "ϕ_max: $ϕ_max"
    
    
    no_P, no_T = length(ES.P), length(ES.T)
    m_k = 0 
    π_avg_k = 0
    k_start = 1
    objective_values = []

    if ~isempty(ES.warmstart)
        warmstart!(ES, ES.warmstart[:method])
        m_k = ES.warmstart[:values][:m]
        π_avg_k = ES.warmstart[:values][:π_avg]
        k_start = 2
    end

    model = Model()
    build_model(model, ES, Dict(Symbol("PDE_type") => PDE_type), model_type="SCP")
    m, π_avg, γ = model[:m], model[:π_avg], model[:γ]
    set_solver(model, "SCP", config_dict_solver)

    global model_dict = Dict()
    global sol_dict = Dict()
    # global violation_dict = Dict()

    # Registered expressions
    global obj_func = objective_function(model)

    if order == 2
        @variable(model, τ[p in ES.P, t in ES.T] >= 0) # auxiliary variable 
    end

    δ = [zeros(no_P,no_T), δ_init/Cost_base*ones(no_P,no_T)]
    γ_k = 0

    # δ_γ_init = 1/(Cost_base*Γ_base)
    δ_γ_init = 1

    start = Dates.now()
    for k in k_start:k_max

        @info "Start iteration $k."
        if k == 1
            # Penalize error term for first iteration
            set_objective_function(model, obj_func+δ_γ_init*sum(γ[p,t]^2 for p in ES.P, t in ES.T))
        else
            # Add Taylor-series-based constraint for momentum conservation constraint
            add_taylor_series_expansion!(model, k, ES, m_k[:,1:end], π_avg_k[:,1:end], order)
            # Add augmentation terms for iteration k
            # set_objective_function(model, 
            #     obj_func + sum(objective_augmentation.(δ[k], π_avg, π_avg_k, m, m_k)))

            if augmentation_method == "L2"
                set_objective_function(model, obj_func + #map_coefficients_inplace!(a -> round(a, digits=8),
                    sum(δ[k][p,t]*((π_avg[p,t] - π_avg_k[p,t])^2 + (m[p,t] - m_k[p,t])^2)
                        for p in ES.P, t in ES.T))
            elseif augmentation_method == "L1"
                if k == 2
                    @variable(model, κ_m[p in ES.P, t in ES.T])
                    @variable(model, κ_π[p in ES.P, t in ES.T])
                    set_objective_function(model, obj_func + sum(δ[k][p,t]*(κ_m[p,t] + κ_π[p,t]) for p in ES.P, t in ES.T))

                    add_L1_norm_constraint(model, ES, δ[k], m_k, π_avg_k)
                else
                    change_L1_norm_constraint_rhs(model, ES, m_k, π_avg_k)
                end
            else throw(ArgumentError("""Passed augmentation method is $(augmentation_method).
                Must be one of 'L1', 'L2'."""))
            end

        end
        # Solve model and extract required optimial values
        solve_model!(model)
        model_dict[k] = copy(model)
        π_avg_sol, γ_sol, m_sol = value.(π_avg), value.(γ), value.(m)
        append!(objective_values, objective_value(model))

        q_d_cur_sol = value.(model[:q_d_cur])
        sol_dict[k] = Dict(:π_avg => π_avg_sol, :γ => γ_sol, :m => m_sol, :q_d_cur => q_d_cur_sol)

        # Check solution convergence with γ from last iteration
        Γ_sgn_max = [m_sol[p,t] == 1 ? ES.pipes[p].Γ_max : ES.pipes[p].Γ_min for p in ES.P, t in ES.T]
        # ϕ = abs.((γ_k .- m_sol[:,2:end].*abs.(m_sol[:,2:end])./π_avg_sol[:,2:end])./Γ_max)
        ϕ = abs.(calc_rel_relaxation_gap.(m_sol[:,1:end], π_avg_sol[:,1:end], γ_sol, Γ_sgn_max))

        # violation_dict[k] = ϕ
        @info "obj:", objective_value(model)
        if all(ϕ .<= ϕ_max)
            write_solution!(ES, model, ES.model_type)
            ES.obj_val = objective_values # overwrite objective value with list including all iterations
            @info "Process converged after $k out of $k_max iterations."
            break
        end

        println("Maximum violation is $(maximum(ϕ)) > $ϕ_max in iteration $k.")

        # Save optimal values for next iteration
        π_avg_k, γ_k, m_k = π_avg_sol, γ_sol, m_sol

        println("End iteration $k.")
        if (k < k_max) && (k > 1)

            # penalization_cost = sum(δ[k][p,t]*(
            #     (π_avg_sol[p,t] - π_avg_k[p,t])^2 + (m_sol[p,t] - m_k[p,t])^2)
            #     for p in ES.P, t in ES.T)
            # y = (penalization_cost+sum(ϕ))/(objective_value(model)+sum(ϕ))

            # convergence_rate_k = min(convergance_rate(y, 1, 1), 10)
            # # update = penalizer_update.(ϕ, 2, 2) * convergence_rate_k/10 .* δ[k]

            update = δ_update*δ[k]
            push!(δ, min.(update, ones(no_P, no_T)*δ_max))
        elseif k == k_max
            println("Process did NOT converge after $k_max iterations.
            The maximum violation is $(maximum(ϕ)) > $ϕ_max.
            You may want to try a higher value of k_max.")
            # include possible metric to capture gap here.
        end

    end
    stop = Dates.now()
    ES.solve_time = ((stop-start)/Millisecond(1000))
    @info "Solution time: $(ES.solve_time) seconds."

    return nothing
end


function delete_and_unregister!(model, cons_names)
    for con_name in cons_names
        delete.(model, model[con_name])
        unregister.(model, con_name)
    end
    return nothing
end


function add_L1_norm_constraint(model, ES, δ, m_k, π_avg_k)
    m, π_avg, κ_m, κ_π = model[:m], model[:π_avg], model[:κ_m], model[:κ_π]

    @constraint(model, L1_m_pos[p in ES.P, t in ES.T], -κ_m[p,t] + m[p,t] <= m_k[p,t])
    @constraint(model, L1_m_neg[p in ES.P, t in ES.T], κ_m[p,t] + m[p,t] >= m_k[p,t])

    @constraint(model, L1_π_pos[p in ES.P, t in ES.T], -κ_π[p,t] + π_avg[p,t] <= π_avg_k[p,t])
    @constraint(model, L1_π_neg[p in ES.P, t in ES.T], κ_π[p,t] + π_avg[p,t] >= π_avg_k[p,t])

    return nothing
end

function change_L1_norm_constraint_rhs(model, ES, m_k, π_avg_k)

    set_normalized_rhs.(model[:L1_m_pos], m_k)
    set_normalized_rhs.(model[:L1_m_neg], m_k)
    set_normalized_rhs.(model[:L1_π_pos], π_avg_k)
    set_normalized_rhs.(model[:L1_π_pos], π_avg_k)

    return nothing
end


function objective_augmentation(δ, π_avg, π_avg_k, m, m_k)
    return δ*((π_avg - π_avg_k)^2 + (m - m_k)^2)
end


function add_taylor_series_expansion!(model, k, ES, m_k, π_avg_k, order)
    # if order == 1
    #     add_first_order_taylor_polynomial!(model, ES, m_k, π_avg_k)
    # else
    #     add_second_order_taylor_polynomial!(model, ES, m_k, π_avg_k)
    # end

    if k > 2
        modify_first_order_taylor_polynomial!(model, m_k, π_avg_k)
    else
        add_first_order_taylor_polynomial!(model, ES, m_k, π_avg_k)
    end
    return nothing

end


function add_first_order_taylor_polynomial!(model, ES, m_k, π_avg_k)
    """
    Add first-order taylor series approximation for the non-convex term in the
    conservation of momentum constraint γ = m*|m|/π_avg.
    In iteration k, this is
    γ = m^k*|m^k|/π_avg^k + 2*|m^k|/π_avg^k*(m-m^k) - m^k*|m^k|/(π_avg^k)^2*(π_avg-π_avg^k),
    which can reformulated as: γ = 2*|m^k|/π_avg^k*(m) - m^k*|m^k|/(π_avg^k)^2*(π_avg)
    γ
    """
    m, π_avg, γ = model[:m], model[:π_avg], model[:γ]
    @constraint(model, taylor_series_expansion[p in ES.P, t in ES.T], γ[p,t] +
        coeff_m(m_k[p,t], π_avg_k[p,t]) * m[p,t] +
        coeff_π_avg(m_k[p,t], π_avg_k[p,t]) * π_avg[p,t] == 0)
    return nothing

end


function coeff_m(m_k, π_avg_k)
    return -2*abs(m_k)/π_avg_k
end


function coeff_π_avg(m_k, π_avg_k)
    return m_k*abs(m_k)/π_avg_k^2
end


function modify_first_order_taylor_polynomial!(model, m_k, π_avg_k)
    m, π_avg, taylor_series_expansion = model[:m], model[:π_avg], model[:taylor_series_expansion]
    set_normalized_coefficient.(
        taylor_series_expansion, m[:,1:end], coeff_m.(m_k, π_avg_k))
    set_normalized_coefficient.(
        taylor_series_expansion, π_avg[:,1:end], coeff_π_avg.(m_k, π_avg_k))
    return nothing
end


function calculate_hessian_γ(m_k, π_avg_k)
    return reshape(
        [2*sign(m_k)/π_avg_k,
         -2*abs(m_k)/π_avg_k^2,
         -2*abs(m_k)/π_avg_k^2,
         2*m_k*abs(m_k)/π_avg_k^3],
            2, 2)
end


function calculate_psd_part_of_matrix(H)
    """
    Calculates the positive semi-definite part of a matrix.
    Assume 
    """
    eig_vals, eig_vecs = eigen(H)
    eig_vals[eig_vals .< 0] .= 0
    # eig_vals = round.(eig_vals, digits=6)
    return eig_vecs * diagm(eig_vals) * eig_vecs'
end


function calculate_psd_hessian(m_k, π_avg_k)
    H = calculate_hessian_γ(m_k, π_avg_k)
    return calculate_psd_part_of_matrix(H)
end


function tsp_sop(m, m_k, π_avg, π_avg_k)
    return 1/2 * [m - m_k, π_avg - π_avg_k]' * 
        calculate_psd_hessian(m_k, π_avg_k) * 
        [m - m_k, π_avg - π_avg_k]
end


function tsp_fop(m, m_k, π_avg, π_avg_k)
    return (abs(m_k) * m_k)/π_avg_k +
        2*abs(m_k)/π_avg_k * (m - m_k) -
        m_k*abs(m_k)/π_avg_k^2 * (π_avg - π_avg_k)
end



function add_second_order_taylor_polynomial!(m, ES, m_k, π_avg_k)
    """
    Add description.
    """
    # Add Taylor-series expansion for current iteration

    m, π_avg, γ, τ = model[:m], model[:π_avg], model[:γ], model[:τ]
    @constraint(m, taylor_series_expansion[p in ES.P, t in ES.T], γ[p,t] ==
        # map_coefficients_inplace!(a -> a <= 1e-10 ? 0 : a,
        tsp_fop(m[p,t], m_k[p,t], π_avg[p,t], π_avg_k[p,t]) + τ[p,t])
    @constraint(m, tsp_second_order_polynomial[p in ES.P, t in ES.T],
        # map_coefficients_inplace!(a -> a <= 1e-10 ? 0 : a,
        tsp_sop(m[p,t], m_k[p,t], π_avg[p,t], π_avg_k[p,t]) <= τ[p,t])
    # @constraint(m, tsp_second_order_polynomial[p in ES.P, t in ES.T],
    #     1/2 * [m[:m][p,t] - m_k[p,t], m[:π_avg][p,t] - π_avg_k[p,t]]' * 
    #     calculate_psd_hessian(m_k[p,t], π_avg_k[p,t]) * 
    #     [m[:m][p,t] - m_k[p,t], m[:π_avg][p,t] - π_avg_k[p,t]] <= m[:τ][p,t]
    # ) 

    return 
end


function penalizer_update(method)
    if method == "constant"
        return constant_update(a)
    elseif method == "log10"
        return log_update(x, a, b)
    end
end


function constant_update(a)
    return a
end


function log_update(x, a, b)
    return -a*log10(x)+b
end


function convergance_rate(y, c, d)
    return -c*log10(y)+d
end