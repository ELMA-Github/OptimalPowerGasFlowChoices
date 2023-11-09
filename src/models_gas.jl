
#=
Which models exist:
PDE
SLP (TSP for γ)
MISOCP (relaxation of γ)
MIQP (gurobi sbnb for γ)
MILP (squared pressures)
MISOCP (sbnb tightening)
PWL (MILP)
=#

# function GS_convex_PDE_model(
#     model::JuMP.Model,
#     nodes::Vector{Node},
#     I::Vector{Int},
#     sources::Vector{GasSource},
#     S::Vector{Int},
#     demands::Vector{GasLoad},
#     D::Vector{Int},
#     pipes::Vector{Pipeline},
#     P::Vector{Int},
#     compressors::Vector{Compressor},
#     C::Vector{Int},
#     dt::Float,
#     T::Vector{Int};
#     PDE_type::Int=3)
#     """
#     Variables:
#     q_s: Supply of gas of source s at timestep t [kg/s].
#     q_d_cur: Curtailment of gas demand d at timestep t [kg/s].
#     π: Pressure at node i at timestep t [MPa].
#     π_avg: Average pressure in pipeline p at timestep t [MPa].
#     m: Average gas flow in pipeline p at timestep t [kg/s].
#     m_in: Gas flow entering pipeline p at timestep t [kg/s].
#     m_out: Gas flow exiting pipeline p at timestep t [kg/s].
#     m_c: Flow of gas through compressor c at timestep t [kg/s].
#     γ: Auxiliary variable used for convexification of conservation of momentum
#         constraint for pipeline p at timestep t.
#         Based on flow defined on two directions.

#     compressor_limit: Ensures that the compression ratio, i.e.,
#         inlet divided by outlet pressure, is within limits
#         for each compressor c at timestep t.
#     """

#     ut1, ut2 = PDE_type_parameters(PDE_type)

#     @variables model begin
#         sources[s].Q_s_min <= q_s[t,s] <= sources[s].Q_s_max
#         0 <= q_d_cur[t,d] <= demands[d].Q_d[t]
#         nodes[i].Π_min <= π[t,i] <= nodes[i].Π_max
#         π_avg[t,p] # bounds?
#         m[p in P, t in T]
#         m_in[p in P, t in T]
#         m_out[p in P, t in T]
#         γ[p in P, t in T]
#     end

#     if ~isempty(compressors)
#         @variable(model, [t in T, c in C], m_c[t,c] >= 0)
#         @constraints model begin
#             compressor_limit_l[c in C, t in T],
#                 compressors[c].CR_min*π[compressors[c].start,t] <= π[compressors[c].stop,t] 
#             compressor_limit_u[c in C, t in T],
#                 π[compressors[c].stop,t] <= compressors[c].CR_max*π[compressors[c].start,t]
#         end
#     end

#     constraint_linepack_end_condition(model, demands, D, S, T)
#     constraint_nodal_gas_balance(model, I, demands, D, pipes, P, compressors, C, T)
#     constraint_average_flow_and_pressure(model, pipes, P, T)
#     constraint_conservation_mass(model, pipes, P, ut1, dt, T)
#     constraint_conservation_momentum(model, pipes, P, ut2, dt, T)

#     return nothing
# end


function build_gas_model(
    model_type::String,
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    P_hat::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    boundary_conditions::Dict{Symbol,Vector},
    dt::Number,
    T::Vector{Int};
    config_dict_algorithm::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

    if model_type == "NCNLP"
        GS_nonconvex_PDE_model(model, nodes, I, sources, S, demands, D, pipes, P, P_hat,
            compressors, C, boundary_conditions, dt, T;
            config_dict_algorithm...)
    elseif model_type == "SCP"
        GS_SLP_model(model, nodes, I, sources, S, demands, D, pipes, P,
            compressors, C, boundary_conditions, dt, T;
            config_dict_algorithm...)            
    # elseif model_type == "nonconvex_lp_form"
    #     ES.obj_val, ES.solve_time, ES.Q_w, ES.Q_nm, ES.Q_nm_in, ES.Q_nm_out,
    #         ES.Qd_cur_gas, ES.Q_c, ES.pr, ES.LP = nonconvex_linepack(ES)
    elseif model_type == "MINLP"
        GS_MINLP_PDE_model(model, nodes, I, sources, S, demands, D, pipes, P, P_hat,
            compressors, C, boundary_conditions, dt, T;
            config_dict_algorithm...)
    elseif model_type == "PWL"
        GS_PWL_PDE_model(model, nodes, I, sources, S, demands, D, pipes, P,
            compressors, C, boundary_conditions, dt, T;
            config_dict_algorithm...)
    elseif model_type in ["MISOCP", "MISOCP_sbnb"]
        GS_MISOCP_PDE_model(model, nodes, I, sources, S, demands, D, pipes, P, P_hat,
            compressors, C, boundary_conditions, dt, T;
            config_dict_algorithm...)
    elseif model_type == "MILP"
        GS_MILP_PDE_model(model, nodes, I, sources, S, demands, D, pipes, P, P_hat,
            compressors, C, boundary_conditions, dt, T;
            config_dict_algorithm...)
    elseif model_type == "MISOCP_McCormick"
        GS_MISOCP_PDE_McCormick_model(model, nodes, I, sources, S, demands, D, pipes, P, P_hat,
            compressors, C, boundary_conditions, dt, T;
            config_dict_algorithm...)
    elseif model_type == "MISOCP_weymouth"
        GS_MISOCP_weymouth(model, nodes, I, sources, S, demands, D, pipes, P, P_hat,
            compressors, C, dt, T;
            config_dict_algorithm...)
    elseif model_type == "CELP"
        GS_CE_PDE_model(model, nodes, I, sources, S, demands, D, pipes, P,
        compressors, C, boundary_conditions, dt, T;
            config_dict_algorithm...)
    elseif model_type == "trade"
        GS_trade_model(model, I, sources, S, demands, D, pipes, P, compressors, C, T)
    elseif model_type == "boundary_conditions"
        GS_boundary_PDE_model(model, nodes, I, sources, S, demands, D, pipes, P, P_hat,
        compressors, C, dt, T;
            config_dict_algorithm...)
    else 
        throw(ArgumentError("""Passed model_type is $(ES.model_type).
        Must be one of 'NCNLP', 'SCP', 'MINLP', 'MISOCP', 'MISOCP_sbnb', 'PWL', 'trade'."""))
    end
    return nothing
end


function GS_boundary_PDE_model(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    P_hat::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    dt::Number,
    T::Vector{Int};
    PDE_type::Int=3)

    ut1, ut2 = PDE_type_parameters(PDE_type)

    variable_gas_supply(model, sources, S, T)
    variable_gas_curtailment(model, demands, D, T)
    variable_nodal_pressure(model, nodes, I, T)
    variable_pipeline_pressure(model, pipes, P, T)
    variable_pipeline_flow_twoway(model, pipes, P, T)
    variable_convex_momentum_auxiliary(model, P, T)

    if ~isempty(compressors)
        variable_compressor_flow(model, C, T)
        constraint_compressor_limits(model, compressors, C, T)
        constraint_compressor_fuel_consumption(model, compressors, C, T)
    end

    expression_nodal_gas_balance(model, I, sources, S, demands, D, pipes, P, compressors, C, T)
    constraint_average_flow_and_pressure(model, pipes, P, T)
    constraint_conservation_mass(model, pipes, P, ut1, dt, T)
    constraint_conservation_momentum(model, pipes, P, ut2, dt, T)

    constraint_linepack_end_condition(model, pipes, P, T)
    constraint_pressure_slack(model, nodes, T)

    constraint_conservation_momentum_nonconvex_twoway(model, P, T)
    constraint_steady_state_initial(model, P)

    return nothing
end


function GS_convex_PDE_model(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    boundary_conditions::Dict{Symbol,Vector},
    dt::Number,
    T::Vector{Int};
    PDE_type::Int=3)

    ut1, ut2 = PDE_type_parameters(PDE_type)

    variable_gas_supply(model, sources, S, T)
    variable_gas_curtailment(model, demands, D, T)
    variable_nodal_pressure(model, nodes, I, T)
    variable_pipeline_pressure(model, pipes, P, T)
    variable_pipeline_flow_twoway(model, pipes, P, T)
    variable_convex_momentum_auxiliary(model, P, T)

    if ~isempty(compressors)
        variable_compressor_flow(model, C, T)
        constraint_compressor_limits(model, compressors, C, T)
        constraint_compressor_fuel_consumption(model, compressors, C, T)
    end

    expression_nodal_gas_balance(model, I, sources, S, demands, D, pipes, P, compressors, C, T)
    constraint_average_flow_and_pressure(model, pipes, P, T)
    constraint_conservation_mass(model, pipes, P, ut1, dt, T)
    constraint_conservation_momentum(model, pipes, P, ut2, dt, T)
    if PDE_type != 1
        constraint_linepack_end_condition_boundary(model, pipes, P, T)
    end
    constraint_pressure_slack(model, nodes, T)

    if PDE_type != 1
        constraint_boundary_conditions(model, boundary_conditions, P)
    end

    return nothing
end


function GS_nonconvex_PDE_model(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    P_hat::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    boundary_conditions::Dict{Symbol,Vector},
    dt::Number,
    T::Vector{Int};
    PDE_type::Int=3)

    GS_convex_PDE_model(
        model, nodes, I, sources, S, demands, D, pipes, P,
        compressors, C, boundary_conditions, dt, T;
        PDE_type)
        
    # constraint_subpipeline_flow(model, pipes, P_hat, T)
    constraint_conservation_momentum_nonconvex_twoway(model, P, T)
    return nothing
end


function GS_SLP_model(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    boundary_conditions::Dict{Symbol,Vector},
    dt::Number,
    T::Vector{Int};
    PDE_type::Int=3)

    GS_convex_PDE_model(
        model, nodes, I, sources, S, demands, D, pipes, P,
        compressors, C, boundary_conditions, dt, T;
        PDE_type)

    return nothing
end


function GS_MINLP_PDE_model(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    P_hat::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    boundary_conditions::Dict{Symbol,Vector},
    dt::Number,
    T::Vector{Int};
    PDE_type::Int=3,
    M::Number=1000
    )

    GS_convex_PDE_model(
        model, nodes, I, sources, S, demands, D, pipes, P,
        compressors, C, boundary_conditions, dt, T;
        PDE_type)

    variable_pipeline_flow_oneway(model, pipes, P, P_hat, T)
    variable_convex_momentum_auxiliary_oneway(model, pipes, P, T) 
    constraint_oneway_flow_limits(model, pipes, P, T)
    constraint_oneway_conservation_momentum(model, pipes, P, T)
    constraint_conservation_momentum_nonconvex_oneway(model, P, T) 

    return nothing
end


function GS_MISOCP_PDE_model(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    P_hat::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    boundary_conditions::Dict{Symbol,Vector},
    dt::Number,
    T::Vector{Int};
    PDE_type::Int=3,
    linear_overestimator::Bool=true)

    GS_convex_PDE_model(
        model, nodes, I, sources, S, demands, D, pipes, P,
        compressors, C, boundary_conditions, dt, T;
        PDE_type)

    variable_pipeline_flow_oneway(model, pipes, P, P_hat, T)
    variable_convex_momentum_auxiliary_oneway(model, pipes, P, T) 
    constraint_oneway_flow_limits(model, pipes, P, T)
    constraint_oneway_conservation_momentum(model, pipes, P, T)
    if linear_overestimator
        constraint_linear_overestimator(model, pipes, P, nodes, T)
    end

    constraint_conservation_momentum_soc(model, P, T)

    return nothing
end


function GS_MILP_PDE_model(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    P_hat::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    boundary_conditions::Dict{Symbol,Vector},
    dt::Number,
    T::Vector{Int};
    PDE_type::Int=3,
    linear_overestimator::Bool=true)

    GS_convex_PDE_model(
        model, nodes, I, sources, S, demands, D, pipes, P,
        compressors, C, boundary_conditions, dt, T;
        PDE_type)

    variable_pipeline_flow_oneway(model, pipes, P, P_hat, T)
    variable_convex_momentum_auxiliary_oneway(model, pipes, P, T) 
    constraint_oneway_flow_limits(model, pipes, P, T)
    constraint_oneway_conservation_momentum(model, pipes, P, T)
    if linear_overestimator
        constraint_linear_overestimator(model, pipes, P, nodes, T)
    end

    constraint_linear_relaxation_conservation_momentum(model, pipes, P, nodes, T)

    return nothing
end


function GS_MISOCP_PDE_McCormick_model(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    P_hat::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    boundary_conditions::Dict{Symbol,Vector},
    dt::Number,
    T::Vector{Int};
    PDE_type::Int=3,
    M::Number=1000)

    GS_convex_PDE_model(
        model, nodes, I, sources, S, demands, D, pipes, P,
        compressors, C, boundary_conditions, dt, T;
        PDE_type)

    variable_pipeline_flow_oneway(model, pipes, P, P_hat, T)
    variable_convex_momentum_auxiliary_oneway(model, pipes, P, T) 
    constraint_oneway_flow_limits(model, pipes, P, T)
    constraint_oneway_conservation_momentum(model, pipes, P, T)
    # constraint_linear_overestimator(model, pipes, P, nodes, T)

    variable_McCormick_oneway(model, P, T) 
    constraint_conservation_momentum_soc_McCormick(model, pipes, P, T)

    # constraint_conservation_momentum_soc(model, P, T)


    return nothing
end


function GS_MISOCP_weymouth(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    P_hat::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    dt::Number,
    T::Vector{Int})
    """
    Model based on Schwele, A., et al. (2020). Coordination of Power and Natural Gas Systems:
    Convexification Approaches for Linepack Modeling,
    Proceedings of IEEE PES PowerTech 2019, DOI: 10.1109/PTC.2019.8810632
    """

    variable_gas_supply(model, sources, S, T)
    variable_gas_curtailment(model, demands, D, T)

    variable_pipeline_flow_twoway(model, pipes, P, T)
    variable_pipeline_flow_oneway(model, pipes, P, P_hat, T)
    # variable_pipeline_flow_oneway_AS(model, pipes, P, T)


    # Special variables for reformulated model
    variable_squared_average_pressure_auxilaries_AS(model, nodes, pipes, P, T) 
    variable_linepack(model, P, T) 

    if ~isempty(compressors)
        variable_compressor_flow_AS(model, compressors, C, nodes, T)
        constraint_compressor_limits_AS(model, compressors, C, T)
        constraint_compressor_fuel_consumption(model, compressors, C, T)
    end

    expression_nodal_gas_balance(model, I, sources, S, demands, D, pipes, P, compressors, C, T)
    constraint_average_flow_AS(model, P, T)
    constraint_linepack_AS(model, pipes, P, dt, T)
    # constraint_average_flow_and_pressure(model, pipes, P, T)
    # constraint_conservation_mass(model, pipes, P, ut1, dt, T)
    # constraint_conservation_momentum(model, pipes, P, ut2, dt, T)
    constraint_McCormick_envelops_AS(model, nodes, pipes, P, T)
    constraint_weymouth_AS(model, pipes, P, T)
    constraint_oneway_flow_limits(model, pipes, P, T)
    # constraint_linear_overestimator(model, pipes, P, T)
    constraint_pressure_slack(model, nodes, T)

    return nothing
end


function GS_trade_model(
    model::JuMP.Model,
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    T::Vector{Int})

    variable_gas_supply(model, sources, S, T)
    variable_gas_curtailment(model, demands, D, T)
    variable_pipeline_flow_without_linepack(model, P, T)
    if ~isempty(compressors)
        variable_compressor_flow(model, C, T)
    end
    constraint_nodal_gas_trade_balance(
        model, I, sources, S, demands, D, pipes, P, compressors, C, T)

    return nothing
end


function GS_CE_PDE_model(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    boundary_conditions::Dict{Symbol,Vector},
    dt::Real,
    T::Vector{Int};
    PDE_type::Int=3)
    """
    CE: Convex envelop
    """

    GS_convex_PDE_model(
        model, nodes, I, sources, S, demands, D, pipes, P,
        compressors, C, boundary_conditions, dt, T;
        PDE_type)
    # variable_convex_envelops_momentum_conservation_auxiliaries(model, P, T)
    constraint_convex_envelops_conservation_momentum(model, pipes, P, nodes, T)
    return nothing
end

    
function GS_PWL_PDE_model(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    I::Vector{Int},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    demands::Dict{Int, GasLoad},
    D::Vector{Int},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    boundary_conditions::Dict{Symbol,Vector},
    dt::Number,
    T::Vector{Int};
    PDE_type::Int=3,
    method::String="ldcc",
    no_m_set::Int=3,
    no_π_set::Int=3)

    """
    Based on Geißler, B. et al. (2012). Using piecewise linear functions for solving MINLPs.    
    """

    @assert method in ["ldcc"] 
    @assert no_m_set >= 3
    @assert no_π_set >= 3

    m_set = [collect(LinRange(pipes[p].M_min, pipes[p].M_max, no_m_set))
        for p in P, t in T]
    π_avg_set = [collect(LinRange(pipes[p].Π_avg_min, pipes[p].Π_avg_max, no_π_set))
        for p in P, t in T]

    # Set number of vertices to 3 as we are using a 2-simplex
    no_vertices = 3

    function recover_γ((m,π)::Tuple{Float64,Float64})
        return m*abs(m)/π
    end

    function combine_grid_points((m,π)::Tuple{Float64,Float64}, γ::Float64)
        return tuple(m,π,γ)
    end

    no_simplices = 2*(no_m_set-1)*(no_π_set-1) # Number of simplices per pipeline and timestep
    no_bits = Int(ceil(log2(no_simplices))) # Number of bits needed to represent simplices
    # Simplices ids in bit representation
    simplices_ids = [bitstring(i)[end-no_bits+1:end] for i in 0:no_simplices-1]

    c = 1 # simplex number
    simplices_mat = Array{Dict}(undef, length(P), length(T)) # Simplices container
    for p in P, t in T
        # Calculate grid points for each pipeline and timestep
        m_π_grid = collect(Base.product(m_set[p,t], π_avg_set[p,t]))
        grid = combine_grid_points.(m_π_grid, recover_γ.(m_π_grid))

        # Create simplices for each pipeline and timestep
        simplices = Dict{String, Vector{Tuple{Float64,Float64,Float64}}}()
        for (i,j) in Base.product(1:no_m_set-1, 1:no_π_set-1)
            coords_vertices = grid[i:i+1, j:j+1]
            simplices[bitstring((c-1)*2)[end-no_bits+1:end]] = coords_vertices[1:end-1]
            simplices[bitstring((c-1)*2+1)[end-no_bits+1:end]] = coords_vertices[2:end]
            c += 1
        end
        simplices_mat[p,t] = simplices
    end

    # Optimization model
    GS_convex_PDE_model(
        model, nodes, I, sources, S, demands, D, pipes, P,
        compressors, C, boundary_conditions, dt, T;
        PDE_type)
    variable_PWL_ldcc(model, no_bits, no_vertices, simplices_ids, P, T)
    constraint_PWL_ldcc(model, no_bits, no_vertices, simplices_mat, simplices_ids, P, T)

    return nothing
end

      
# function GS_steady_state_model(
#     model::JuMP.Model,
#     nodes::Vector{Node},
#     I::Vector{Int},
#     sources::Vector{GasSource},
#     S::Vector{Int},
#     demands::Vector{GasLoad},
#     D::Vector{Int},
#     pipes::Vector{Pipeline},
#     P::Vector{Int},
#     compressors::Vector{Compressor},
#     C::Vector{Int},
#     dt::Float,
#     T::Vector{Int})
#     """
#     Does not consider linepack nor transients (ut1=0, ut2=0). Therefore, pressures can be substituted to
#     squared pressures, which can then be further formulated.
#     """
#     variable_gas_supply(model, sources, S, T)
#     variable_gas_curtailment(model, demands, D, T)
#     variable_nodal_pressure_squared(model, nodes, I, T)
#     variable_pipeline_pressure(model, pipes, P, T)
#     variable_pipeline_flow_twoway(model, P, T)
#     variable_convex_momentum_auxiliary(model, P, T)

#     if ~isempty(compressors)
#         variable_compressor_flow(model, C, T)
#         constraint_compressor_limits(model, compressors, C, T)
#     end

#     constraint_linepack_end_condition(model, demands, D, S, T)
#     constraint_nodal_gas_balance(model, I, demands, D, pipes, P, compressors, C, T)
#     constraint_average_flow_and_pressure(model, pipes, P, T)
#     constraint_conservation_mass(model, pipes, P, ut1, dt, T)
#     constraint_conservation_momentum(model, pipes, P, ut2, dt, T)

#     return nothing
# end