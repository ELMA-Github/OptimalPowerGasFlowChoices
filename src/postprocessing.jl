include("objective.jl")
################### Calculating pressures ####################
# TODO:Rework this function!
function calc_pressures_from_ϕ_p(
    ϕ⁺::Matrix{<:Number},#::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    ϕ⁻::Matrix{<:Number},#::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    I::Vector{Int})
    # π = Array{<:Number}(undef, maximum(I), size(ϕ⁺)[2])
    π = zeros(maximum(I), size(ϕ⁺)[2])
    for p in P
        π[pipes[p].start,:] = 1/2*(ϕ⁺[p,:].+ϕ⁻[p,:])
        π[pipes[p].stop,:] = 1/2*(ϕ⁺[p,:].-ϕ⁻[p,:])
    end
    return π
end


function calc_pressures_from_ϕ_c(
    ϕ⁺_c::Matrix{<:Number},#::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    ϕ⁻_c::Matrix{<:Number},#::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    compressors::Dict{Int,Compressor},
    C::Vector{Int})
    π = Dict{Int,Vector{<:Number}}()
    for c in C
        π[compressors[c].start] = 1/2*(ϕ⁺_c[c,:].+ϕ⁻_c[c,:])
        π[compressors[c].stop] = 1/2*(ϕ⁺_c[c,:].-ϕ⁻_c[c,:])
    end
    return π
end


function calc_pressurs_from_ϕ(
    ϕ⁺::Matrix{<:Number},
    ϕ⁻::Matrix{<:Number},
    ϕ⁺_c::Matrix{<:Number},
    ϕ⁻_c::Matrix{<:Number},
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    I::Vector{Int},
    compressors::Dict{Int,Compressor},
    C::Vector{Int}
)
    π = calc_pressures_from_ϕ_p(ϕ⁺, ϕ⁻, pipes, P, I)
    π_c = calc_pressures_from_ϕ_c(ϕ⁺_c, ϕ⁻_c, compressors, C)
    for c in C
        π[compressors[c].start,:] = π_c[compressors[c].start]
        π[compressors[c].stop,:] = π_c[compressors[c].stop]
    end
    return π
end


#################### Relaxation Gap ####################
function calc_abs_relaxation_gap(m::Number, π_avg::Number, γ::Number)
    """
    Calculate the absolute relaxation gap for the non-convex part of the
    conservation of momentum constraint.
    """
    return γ .- m*abs(m)/π_avg 
end


function calc_abs_relaxation_gap(ES::EnergySystem)
    return calc_abs_relaxation_gap.(ES.m, ES.π_avg, ES.γ)
end


function calc_rel_relaxation_gap(m::Number, π_avg::Number, γ::Number, Γ_sgn_max::Number)
    return calc_abs_relaxation_gap(m, π_avg, γ)/Γ_sgn_max 
end


function calc_rel_relaxation_gap(ES::EnergySystem)
    Γ_sgn_max = [ES.z[p,t] == 1 ? ES.pipes[p].Γ_max : ES.pipes[p].Γ_min
        for p in ES.P, t in ES.T]
    return calc_rel_relaxation_gap.(ES.m, ES.π_avg[:,2:end], ES.γ, Γ_sgn_max)
end


#################### Relaxation Gap Metrics ####################
function inf_norm(x::Matrix, P::Vector{Int}, T::Union{Vector{Int}, UnitRange})
    return maximum(abs.(x[P,T]))
end


function inf_norm(x::Array)
    return maximum(abs.(x))
end


function root_mean_squared_gap(x::Matrix, P::Vector{Int}, T::Union{Vector{Int}, UnitRange})
    return sqrt(mean(x[P,T].^2))
end


#################### Linepack Metrics ####################
function cycles(h::Matrix, h_max::Vector, P::Vector{Int}, T::Union{Vector{Int}, UnitRange})
    return mean(sum(abs.((h[P,2:end]-h[P,1:end-1])[:,T])./(2*h_max),dims=2))
end


function total_net_LP_change(h::Matrix, P::Vector{Int}, T::Union{Vector{Int}, UnitRange})
    return sum(abs.((h[P,2:end]-h[P,1:end-1])[:,T]))
end


function root_mean_squared_change_of_charge(
    h::Matrix,
    h_max::Vector,
    P::Vector{Int},
    T::Union{Vector{Int}, UnitRange},
    T_sub::Int,
    collapse_subtimesteps::Bool=true)

    rel_lp_diff = (h[P,2:end]-h[P,1:end-1])[:,T]./h_max
    if collapse_subtimesteps
        rel_lp_diff = reduce(hcat, [sum(
            rel_lp_diff[:,(t-1)*T_sub+1:t*T_sub], dims=2) for t in 1:Int(length(T)/T_sub)])
    end

    return sqrt(mean((rel_lp_diff).^2))
end


function inf_norm_change_of_charge(
    h::Matrix,
    h_max::Vector,
    P::Vector{Int},
    T::Union{Vector{Int}, UnitRange},
    T_sub::Int,
    collapse_subtimesteps::Bool=true)

    rel_lp_diff = (h[P,2:end]-h[P,1:end-1])[:,T]./h_max
    if collapse_subtimesteps
        rel_lp_diff = reduce(hcat, [sum(
            rel_lp_diff[:,(t-1)*T_sub+1:t*T_sub], dims=2) for t in 1:Int(length(T)/T_sub)])
    end

    return inf_norm(rel_lp_diff)
end


#################### Metrics ####################
function calc_metrics(
    ES::EnergySystem, 
    P::Vector{Int},
    hours::Union{Vector{Int},UnitRange}=[],
    collapse_subtimesteps::Bool=true)

    if isempty(hours) hours = 1:ES.timehorizon end
    T = Int((hours[1]-1)*3600/ES.dt+1):Int(hours[end]*3600/ES.dt)
    LP_max = [ES.pipes[p].Π_avg_max*ES.pipes[p].S for p in P]
    gap = calc_rel_relaxation_gap(ES)
    Φ_inf = inf_norm(gap, P, T)
    Φ_RMSG = root_mean_squared_gap(gap, P, T)
    Ξ_cyc = cycles(ES.LP, LP_max, P, T)
    Ξ = total_net_LP_change(ES.LP, P, T)
    Ξ_RMSCoC = root_mean_squared_change_of_charge(
        ES.LP, LP_max, P, T, Int(3600/ES.dt), collapse_subtimesteps)
    Ξ_inf = inf_norm_change_of_charge(
        ES.LP, LP_max, P, T, Int(3600/ES.dt), collapse_subtimesteps)
    stats = [Φ_inf, Φ_RMSG, Ξ, Ξ_cyc, Ξ_RMSCoC, Ξ_inf]
    return stats
end 


# reduce(vcat, [calc_metrics(ES1, [p]) for p in ES1.P]')
# reduce(vcat, [calc_metrics(ES2, [p]) for p in ES2.P]')

# calc_metrics(ES1, ES1.P)
# calc_metrics(ES2, ES2.P)

# function calc_metrics(ES::EnergySystem)
#     stats = mean.(calc_metrics_all_pipes(ES))
#     stats = vcat(ES.model_type*"_"*string(ES.dt), stats)
#     return stats
# end 

# function calc_metrics(ES::EnergySystem, p::Int)
#     stats = calc_metrics_all_pipes(ES)
#     stats = [s[p] for s in stats]
#     stats = vcat(ES.model_type*"_"*string(ES.dt), stats)
#     return stats
# end 

function make_metrics_df(
    ESs::Dict{String,<:EnergySystem},
    P::Vector{Int},
    model_names::Vector{String},
    hours::Union{Vector{Int}, UnitRange}=[],
    collapse_subtimesteps::Bool=true
    )
    if isempty(hours) hours = 1:ESs[model_names[1]].timehorizon end
    metrics_mat = reduce(hcat, [calc_metrics(ESs[mn], P, hours, collapse_subtimesteps) for mn in model_names])

    T_dict = Dict(mn => Int((hours[1]-1)*3600/ESs[mn].dt+1):Int(hours[end]*3600/ESs[mn].dt) for mn in model_names)
    el_load_shed = round.([sum(ESs[mn].p_d_cur[:,T_dict[mn]])*ESs[mn].bases_dict[:S_base]*ESs[mn].dt/3600
        for mn in model_names], digits=2)
    total_costs = [ESs[mn].total_cost for mn in model_names]
    solve_times = round.([ESs[mn].solve_time for mn in model_names], digits=2)

    # models = metrics_mat[1,:]
    metrics_mat = collect(round.(metrics_mat', digits=4))
    metrics_mat = hcat(model_names,metrics_mat,el_load_shed,total_costs,solve_times)

    col_names = [
        :Model, :inf_norm, :rmsg, :tot_net_LP_change, :cyc, :rmscoc, :inf_coc, :shed, :cost, :solve_time]

    df = DataFrame(metrics_mat, col_names)
    return df
end


#################### Calculating linepack ex-post ####################
function calc_linepack_expost(π_avg::Number, S::Number)
    """
    Calculate ex-post linepack based on the average pipeline pressure and the
    linepack parameter S.
    """
    return π_avg * S
end


function calc_linepack_expost(ES::EnergySystem)
    S_matrix = [ES.pipes[p].S for p in ES.P, t in vcat(0,ES.T)]
    return calc_linepack_expost.(ES.π_avg, S_matrix)
end


#################### Pipeline pressure drop ####################
function pipe_pressure_drop(ES, p)
    start, stop = ES.pipes[p].start, ES.pipes[p].stop
    return ES.π[start,:] .- ES.π[stop,:]
end

# #################### Check violation of flow bounds ####################
# function check_possible_flow_bound_violation(val::Number, bound::Number)
#     """
#     Check if the absolute of a value is larger than the absolute of the
#     respective bound minus a small margin.
#     """
#     return abs(val) > (abs(bound)*0.99)
# end




# function check_possible_flow_bound_violation(ES::EnergySystem)
#     """
#     Check if the ex-post flows violate the bounds approximated from the Weymouth Equation
#     when considering the transient term. This indicates that the optimal flow
#     might create a violation of the technically feasible flow.
#     The transient term is only considered for t > 1 in the optimization model.
#     """
#     m = ES.m

#     C, T_base = SPEED_OF_SOUND/ES.bases_dict[:C_base], ES.bases_dict[:T_base]

#     M_max = [ES.pipes[p].M_max for p in ES.P, t in ES.T]
#     M_min = [ES.pipes[p].M_min for p in ES.P, t in ES.T]

#     transient_term = [2*ES.pipes[p].D*ES.pipes[p].A/(ES.dt/T_base*ES.pipes[p].F*C^2) 
#         for p in ES.P, t in ES.T[2:end]] .* (m[:,2:end] .- m[:,1:end-1])
#     transient_term = hcat(zeros(length(ES.P)), transient_term)

#     for bound in [("maximum", M_max), ("minimum", M_min)]
#         if bound[1] == "maximum"
#             mat = check_possible_flow_bound_violation.(max.(m,0), bound[2] .- transient_term)
#         else 
#             mat = check_possible_flow_bound_violation.(min.(m,0), bound[2] .- transient_term)
#         end

#         indices = findall(x -> x == 1, mat)
#         if ~isempty(indices)
#             for idx in indices
#                 @warn "Violation of transient $(bound[1]) pipeline flow in pipeline $(Tuple(idx)[1]) in timestep $(Tuple(idx)[2])."
#             end
#         end
#     end

#     return nothing
# end


#################### Calculating linepack ex-post ####################
function calc_compression_ratio(π_in::Number, π_out::Number)
    """
    Calculate the ratio between output and input pressure at compressors.
    """
    return π_out/π_in
end


function calc_compression_ratio(ES::EnergySystem)
    compressors = [ES.compressors[c] for c in ES.C]
    nodes = [(c.start, c.stop) for c in compressors]
    return reduce(vcat, [calc_compression_ratio.(ES.π[fr,:], ES.π[to,:]) for (fr,to) in nodes]')
end


function calc_compression_diff(π_in::Number, π_out::Number)
    """
    Calculate the difference between output and input pressure at compressors.
    """
    return π_out-π_in
end


function calc_compression_diff(ES::EnergySystem)
    compressors = [ES.compressors[c] for c in ES.C]
    nodes = [(c.start, c.stop) for c in compressors]
    return reduce(vcat, [calc_compression_diff.(ES.π[fr,:], ES.π[to,:]) for (fr,to) in nodes]')
end
