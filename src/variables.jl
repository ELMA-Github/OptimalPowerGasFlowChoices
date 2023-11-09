# PDE model

#################### Gas system variables ####################
function variable_gas_supply(
    model::JuMP.Model, sources::Dict{Int, GasSource}, S::Vector{Int}, T::Vector{Int})
    """
    Supply of gas of source s at timestep t [kg/s].
    """
    @variable(model, sources[s].Q_s_min <= q_s[s in S, t in T] <= sources[s].Q_s_max)
    return nothing
end


function variable_gas_curtailment(
    model::JuMP.Model, demands::Dict{Int, GasLoad}, D::Vector{Int}, T::Vector{Int})
    """
    Curtailment of gas demand d at timestep t [kg/s].
    """
    @variable(model, 0 <= q_d_cur[d in D, t in T] <= demands[d].Q_d[t])
    return nothing
end


# TODO: what about compressor flow limits?
function variable_compressor_flow(
    model::JuMP.Model, C::Vector{Int}, T::Vector{Int})
    """
    m_c # Flow of gas through compressor c at timestep t [kg/s].
    m_cf # Flow of gas consumed by compressor c at node i at timestep t [kg/s].
    """
    @variable(model, m_c[c in C, t in T] >= 0)
    @variable(model, m_cf[c in C, t in T] >= 0)
    return nothing
end


function variable_nodal_pressure(
    model::JuMP.Model, nodes::Dict{Int, Node}, I::Vector{Int}, T::Vector{Int})
    """
    Pressure at node i at timestep t [MPa].
    """
    @variable(model, nodes[i].Π_min <= π[i in I, t in T] <= nodes[i].Π_max)
    return nothing
end


function variable_pipeline_pressure(
    model::JuMP.Model, pipes::Dict{Int, Pipeline}, P::Vector{Int}, T::Vector{Int})
    """
    Average pressure in pipeline p at timestep t [MPa].
    """
    T_ext = vcat(0, T)
    @variable(model, pipes[p].Π_avg_min <= π_avg[p in P, t in T_ext] <= pipes[p].Π_avg_max)
    return nothing
end


#todo: should we make bounds on flow?
function variable_pipeline_flow_twoway(
    model::JuMP.Model, pipes::Dict{Int,Pipeline}, P::Vector{Int}, T::Vector{Int})
    """
    Mass flow m can be positive or negative.
    m # Average gas flow in pipeline p at timestep t [kg/s].
    m_in # Gas flow entering pipeline p at timestep t [kg/s].
    m_out # Gas flow exiting pipeline p at timestep t [kg/s].
    """
    T_ext = vcat(0, T)
    @variables model begin
        pipes[p].M_min <= m[p in P, t in T_ext] <= pipes[p].M_max
        m_in[p in P, t in T]
        m_out[p in P, t in T]
    end
    return nothing
end


function variable_pipeline_flow_without_linepack(
    model::JuMP.Model, P::Vector{Int}, T::Vector{Int})
    """
    Mass flow m can be positive or negative.
    m: Gas flow in pipeline p at timestep t [kg/s].
    """
    @variables model begin
        m[p in P, t in T]
    end
    return nothing
end


function variable_pipeline_flow_oneway_AS(
    model::JuMP.Model, pipes::Dict{Int,Pipeline}, P::Vector{Int}, T::Vector{Int}) 
    """
    Mass flow m⁺ and m⁻ are always non-negative. A non-negative flow in positive direction indicates a flow in
    the same direction as the pipe is defined. For example, a pipe that starts at node 1 and stops in node 2,
    the gas flow m⁺ is non-negative when going from node 1 to 2.
    m⁺ # Gas flow in positive direction pipeline p at timestep t [kg/s].
    m⁺_in # Gas flow in positive direction entering pipeline p at timestep t [kg/s].
    m⁺_out # Gas flow in positive direction exiting pipeline p at timestep t [kg/s].
    m⁻ # Gas flow in positive direction pipeline p at timestep t [kg/s].
    m⁻_in # Gas flow in positive direction entering pipeline p at timestep t [kg/s].
    m⁻_out # Gas flow in positive direction exiting pipeline p at timestep t [kg/s].
    u # Indidicates if flow is in positive or negative direction in pipeline p at timestep t.
        For pipeline p from node i to node j, u=1 if flow is going from i to j and vice versa.
    """
    @variables model begin
        pipes[p].M_min <= m[p in P, t in T] <= pipes[p].M_max
        0 <= m⁺[p in P, t in T] <= pipes[p].M_max
        0 <= m⁺_in[p in P, t in T]
        0 <= m⁺_out[p in P, t in T]
        0 <= m⁻[p in P, t in T] <= -pipes[p].M_min
        0 <= m⁻_in[p in P, t in T] 
        0 <= m⁻_out[p in P, t in T]
        z[p in P, t in T], Bin
    end
    return nothing
end


function variable_squared_average_pressure_auxilaries_AS(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    pipes::Dict{Int,Pipeline},
    P::Vector{Int}, 
    T::Vector{Int}) 
    """
    ϕ⁺ # auxiliary variable for pressures ϕ⁺ = π_m + π_n
    ϕ⁻ # auxiliary variabls for pressures ϕ⁻ = π_m - π_n
    ψ # auxilary variable used for McCormick envelops, ψ = ϕ⁺ * ϕ⁻
    """
    @variables model begin
        ϕ⁺[p in P, t in T]
        ϕ⁻[p in P, t in T]
        ψ[p in P, t in T]
    end
    set_lower_bound.(ϕ⁺, [nodes[pipes[p].start].Π_min + nodes[pipes[p].stop].Π_min for p in P])
    set_upper_bound.(ϕ⁺, [nodes[pipes[p].start].Π_max + nodes[pipes[p].stop].Π_max for p in P])
    set_lower_bound.(ϕ⁻, [nodes[pipes[p].start].Π_min - nodes[pipes[p].stop].Π_max for p in P])
    set_upper_bound.(ϕ⁻, [nodes[pipes[p].start].Π_max - nodes[pipes[p].stop].Π_min for p in P])
    return nothing
end


# TODO: what about compressor flow limits?
function variable_compressor_flow_AS(
    model::JuMP.Model, compressors::Dict{Int,Compressor}, C::Vector{Int}, nodes::Dict{Int, Node}, T::Vector{Int})
    """
    ϕ_c # Auxiliary variable for squared average pressure in compressor c  at timestep t.
    m_c # Flow of gas through compressor c at timestep t [kg/s].
    """
    @variables model begin
        m_c[c in C, t in T] >= 0
        ϕ⁺_c[c in C, t in T]
        ϕ⁻_c[c in C, t in T]
        m_cf[c in C, t in T] >= 0
    end
    set_lower_bound.(ϕ⁺_c, [nodes[compressors[c].start].Π_min + nodes[compressors[c].stop].Π_min for c in C])
    set_upper_bound.(ϕ⁺_c, [nodes[compressors[c].start].Π_max + nodes[compressors[c].stop].Π_max for c in C])
    set_lower_bound.(ϕ⁻_c, [nodes[compressors[c].start].Π_min - nodes[compressors[c].stop].Π_max for c in C])
    set_upper_bound.(ϕ⁻_c, [nodes[compressors[c].start].Π_max - nodes[compressors[c].stop].Π_min for c in C])
    return nothing
end


function variable_linepack(
    model::JuMP.Model, P::Vector{Int}, T::Vector{Int}) 
    """
    LP # Linepack mass
    """
    @variables model begin
        LP[p in P, t in T]
    end
    return nothing
end


function variable_pipeline_flow_oneway(
    model::JuMP.Model, pipes::Dict{Int, Pipeline}, P::Vector{Int}, P_hat::Vector{Int}, T::Vector{Int}) 
    """
    This variables are added onto the existing flow variables m, m_in, and m_out.
    Mass flow m⁺ and m⁻ are always non-negative. A non-negative flow in positive direction indicates a flow in
    the same direction as the pipe is defined. For example, a pipe that starts at node 1 and stops in node 2,
    the gas flow m⁺ is non-negative when going from node 1 to 2.
    m⁺: Gas flow in positive direction pipeline p at timestep t [kg/s].
    m⁻: Gas flow in positive direction pipeline p at timestep t [kg/s].
    u: Indidicates if flow is in positive or negative direction in pipeline p at timestep t.
        For pipeline p from node i to node j, u=1 if flow is going from i to j and vice versa.
    """
    @variables model begin
        0 <= m⁺[p in P, t in T] <= pipes[p].M_max
        0 <= m⁻[p in P, t in T] <= -pipes[p].M_min
        z[p in P, t in T], Bin
    end
    return nothing
end


function variable_convex_momentum_auxiliary_oneway(
    model::JuMP.Model, pipes::Dict{Int, Pipeline}, P::Vector{Int}, T::Vector{Int}) 
    """
    Auxiliary variables used for convexification of the conservation of momentum constraint
    for pipeline p at timestep t. Based on flows defined on a single direction only.
    """
    @variables model begin
        0 <= γ⁺[p in P, t in T] <= pipes[p].Γ_max
        0 <= γ⁻[p in P, t in T] <= -pipes[p].Γ_min
    end
    return nothing
end


function variable_McCormick_oneway(
    model::JuMP.Model, P::Vector{Int}, T::Vector{Int}) 
    """
    Auxiliary variables used for building McCormick envelop of π_avg*γ
    for pipeline p at timestep t. Based on flows defined on a single direction only.
    """
    @variables model begin
        v⁺[p in P, t in T]
        v⁻[p in P, t in T]
    end
    return nothing
end


function variable_convex_momentum_auxiliary(
    model::JuMP.Model, P::Vector{Int}, T::Vector{Int})
    """
    Auxiliary variables used for convexification of the conservation of momentum constraint
    for pipeline p at timestep t. Based on flow defined on two directions.
    """
    @variable(model, γ[p in P, t in T])
    return nothing
end


function variable_PWL_ldcc(
    model::JuMP.Model,
    no_bits::Int,
    no_vertices::Int,
    simplices_ids::Vector{String},
    P::Vector{Int},
    T::Vector{Int})
    """
    Piecewise-linearization based on the logarithmic disaggregated convex combination method from
    Geißler, B. et al. (2012). Using piecewise linear functions for solving MINLPs.
    λ: weights or vertices
    u: indicator for active simplex
    """
    @variables model begin
        λ[p in P, t in T, s_id in simplices_ids, j in 1:no_vertices] >= 0
        u[p in P, t in T, l in 1:no_bits], Bin
    end
    return nothing
end


# function variable_nodal_pressure_squared(
#     model::JuMP.Model, nodes::Dict{Int, Pipeline}, I::Vector{Int}, T::Vector{Int})
#     @variable(model, nodes[i].Π_min^2 <= π_sq[i in I, t in T],  <= nodes[i].Π_max^2)
#     return nothing
# end


# function variable_piece_wise_linearization()
# end

function variable_convex_envelops_momentum_conservation_auxiliaries(
    model::JuMP.Model, P::Vector{Int}, T::Vector{Int})
    """
    Auxiliary variables used for building convex envelops of the nonconvex
    conservation of momentum constraint for pipeline p at timestep t.
    Based on Mhanna, S. et al. (2021). Iterative LP-based Methods for the
    Multiperiod Optimal Electricity and Gas Flow Problem.
    """
    @variables model begin
        u[p in P, t in T]
        v[p in P, t in T]
    end
    return nothing
end

#################### Power system variables ####################
function variable_gas_consumption_gfpp(
    model::JuMP.Model, GFPPs::Vector{Int}, T::Vector{Int})
    """
    Gas consumption of gas-fired power plant g at timestep t [kg/s].
    """
    @variable(model, q_g[g in GFPPs, t in T])
    return nothing
end


function variable_power_curtailment(
    model::JuMP.Model, demands_EL::Dict{Int,PowerDemand}, D_el::Vector{Int}, T::Vector{Int})
    """
    Curtailment of power demand d_el at timestep t [MW].
    """
    @variable(model, 0 <= p_d_cur[d in D_el, t in T] <= demands_EL[d].D_el[t])
    return nothing
end


function variable_power_production_dispatchables(
    model::JuMP.Model, generators::Dict{Int, DispatchableGenerator}, G::Vector{Int}, T::Vector{Int})
    """
    Power production of dispatchable unit g at timestep t [MW].
    """
    @variable(model, 0 <= p[g in G, t in T] <= generators[g].P_max)
    return nothing
end


function variable_power_production_renewables(
    model::JuMP.Model, windgenerators::Dict{Int, WindGenerator}, J::Vector{Int}, T::Vector{Int})
    """
    Power production of renewable unit j at timestep t [MW].
    """
    @variable(model, 0 <= w[j in J, t in T] <= windgenerators[j].P_w[t])
    return nothing
end


function variable_bus_voltage_angles(
    model::JuMP.Model, N::Vector{Int}, T::Vector{Int})
    """
    Voltage phase angle of bus n at timestep t [rad].
    """
    @variable(model, θ[n in N, t in T])
    return nothing
end