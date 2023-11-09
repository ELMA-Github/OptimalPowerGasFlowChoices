###############################################################################
# General constraints common for all models
# constraint_linepack_end_condition
# constraint_compressor_limits
# constraint_nodal_gas_balance

###############################################################################


#################### Gas system constraints ####################
function constraint_boundary_conditions(
    model::JuMP.Model,
    boundary_conditions::Dict{Symbol,Vector},
    P::Vector{Int})
    """
    Fixes the boundary coniditons for the average pressure and flow at the beginning
    of the time horizon to a given value.
    """
    m, π_avg = model[:m], model[:π_avg]
    @constraints model begin
        boundary_condition_flow[p in P],
            m[p,0] == boundary_conditions[:m][p]
        boundary_condition_pressure[p in P],
            π_avg[p,0] == boundary_conditions[:π_avg][p]
    end
    return nothing
end


function constraint_steady_state_initial(
    model::JuMP.Model,
    P::Vector{Int})
    """
    Ensures that average flows and pressures are the same as initially, i.e.,
    in a steady-state.
    """
    m, π_avg = model[:m], model[:π_avg]
    @constraints model begin
        boundary_condition_flow[p in P],
            m[p,0] == m[p,1]
        boundary_condition_pressure[p in P],
            π_avg[p,0] == π_avg[p,1]
    end
    return nothing
end


function constraint_linepack_end_condition(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    T::Vector{Int})
    """
    Ensures that the total amount of gas extracted from the whole gas system is at least also injected.
    Including natural gas-fired power plants in the power system.
    """
    π_avg = model[:π_avg]
    @constraint(model, linepack_end_condition,
        sum(π_avg[p,T[end]]*pipes[p].S for p in P) >= 
            sum(π_avg[p,T[1]]*pipes[p].S for p in P))
    return nothing
end


function constraint_linepack_end_condition_boundary(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    T::Vector{Int})
    """
    Ensures that the total amount of gas extracted from the whole gas system is at least also injected.
    Including natural gas-fired power plants in the power system.
    """
    π_avg = model[:π_avg]
    # @constraint(model, linepack_end_condition,
    #    sum(π_avg[p,T[end]]*pipes[p].S for p in P) >= sum(π_avg[p,0]*pipes[p].S for p in P))
    @constraint(model, linepack_end_condition[p in P], π_avg[p,T[end]] >= π_avg[p,0])
    return nothing
end


function constraint_pressure_slack(
    model::JuMP.Model, nodes::Dict{Int,Node}, T::Vector{Int})
    """
    Fixes pressure at the reference node i_slack to Π_slack at timestep t.
    """
    π = model[:π]
    I_slack = findall(x->isa(x, SlackNode), nodes)
    @constraint(model, pressure_slack_node[i in I_slack, t in T],
        π[i, t] == nodes[i].Π_slack)
    return nothing
end


function constraint_compressor_limits(
    model::JuMP.Model, compressors::Dict{Int,Compressor}, C::Vector{Int}, T::Vector{Int})
    """
    Ensures that the compression ratio, i.e., inlet divided by outlet pressure, is within limits
    for each compressor c at timestep t.
    """
    π = model[:π]
    @constraints model begin
        compressor_limit_l[c in C, t in T],
            compressors[c].CR_min*π[compressors[c].start,t] <= π[compressors[c].stop,t] 
        compressor_limit_u[c in C, t in T],
            π[compressors[c].stop,t] <= compressors[c].CR_max*π[compressors[c].start,t]
    end
    return nothing
end


function constraint_compressor_fuel_consumption(
    model::JuMP.Model, compressors::Dict{Int,Compressor}, C::Vector{Int}, T::Vector{Int})
    """
    Defines the fuel consumption of each compressor c at timestep t as a constant fraction
    of gas flow through the compressor.
    """
    m_c, m_cf, π = model[:m_c], model[:m_cf], model[:π]
    @constraint(model, compressor_fuel_consumption[c in C, t in T],
        m_cf[c,t] == m_c[c,t]*compressors[c].fuel_consumption)# +
            # (π[compressors[c].stop,t]-π[compressors[c].start,t])/54*compressors[c].fuel_consumption/20
        # )
    return nothing
end


#TODO: rework this function
function constraint_compressor_limits_AS(
    model::JuMP.Model,
    compressors::Dict{Int,Compressor},
    C::Vector{Int},
    T::Vector{Int})
    """
    Ensures that the compression ratio, i.e., inlet divided by outlet pressure (represented by
    auxiliary variables), is within limits for each compressor c at timestep t.
    """
    ϕ⁺_c, ϕ⁻_c = model[:ϕ⁺_c], model[:ϕ⁻_c]
    @constraints model begin
        compressor_limit_l[c in C, t in T],
            compressors[c].CR_min*0.5*(ϕ⁺_c[c,t]+ϕ⁻_c[c,t]) <= 0.5*(ϕ⁺_c[c,t]-ϕ⁻_c[c,t])
        compressor_limit_u[c in C, t in T],
            0.5*(ϕ⁺_c[c,t]-ϕ⁻_c[c,t]) <= compressors[c].CR_max*0.5*(ϕ⁺_c[c,t]+ϕ⁻_c[c,t])
    end
    return nothing
end


function constraint_nodal_gas_balance_GS(
    model::JuMP.Model,
    I::Vector{Int},
    T::Vector{Int}
)
    expr_nodal_gas_balance = model[:expr_nodal_gas_balance] 
    @constraint(model, nodal_gas_balance[i in I, t in T], expr_nodal_gas_balance[i,t] == 0)
    return nothing
end


function constraint_nodal_gas_balance_IES(
    model::JuMP.Model,
    I::Vector{Int},
    generators::Dict{Int,DispatchableGenerator},
    G::Vector{Int}, 
    T::Vector{Int}
)
    GFPPs = [g for g in G if isa(generators[g], GasFiredGenerator)]
    expr_nodal_gas_balance, p = model[:expr_nodal_gas_balance], model[:p]
    @constraint(model, nodal_gas_balance[i in I, t in T],
        expr_nodal_gas_balance[i,t] -
        sum(generators[g].Γ*p[g,t] for g in GFPPs if generators[g].node == i) == 0)
    return nothing
end


function expression_nodal_gas_balance(
    model::JuMP.Model,
    I::Vector{Int},
    sources::Dict{Int,GasSource},
    S::Vector{Int},
    demands::Dict{Int,GasLoad},
    D::Vector{Int},
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    compressors::Dict{Int,Compressor},
    C::Vector{Int},
    T::Vector{Int})
    """
    Ensures balance between gas injections and extractions at each node i at time step t.
    Accounts for flows through compressors if there are any.
    """

    q_s, q_d_cur, m_in, m_out = 
        model[:q_s], model[:q_d_cur], model[:m_in], model[:m_out]

    nodal_gas_balance = [
        sum(q_s[s,t] for s in S if sources[s].node == i; init = 0.0) -
        sum(demands[d].Q_d[t] for d in D if demands[d].node == i; init = 0.0) +
        sum(q_d_cur[d,t] for d in D if demands[d].node == i; init = 0.0) +
        sum(m_out[p,t] for p in P if pipes[p].stop == i; init = 0.0) -
        sum(m_in[p,t] for p in P if pipes[p].start == i; init = 0.0)
            for i in I, t in T]

    if isempty(compressors) == true
        @expression(model, expr_nodal_gas_balance[i in I, t in T], nodal_gas_balance[i,t])
    else
        m_c, m_cf = model[:m_c], model[:m_cf]
        @expression(model, expr_nodal_gas_balance[i in I, t in T], 
            nodal_gas_balance[i,t] +
            sum(m_c[c,t] for c in C if compressors[c].stop == i) -
            sum(m_c[c,t] for c in C if compressors[c].start == i) -
            sum(m_cf[c,t] for c in C if compressors[c].fuel_node == i)
            )
    end
    return nothing
end


# function expression_nodal_gas_balance_AS(
#     model::JuMP.Model,
#     I::Vector{Int},
#     sources::Dict{Int,GasSource},
#     S::Vector{Int},
#     demands::Dict{Int,GasLoad},
#     D::Vector{Int},
#     pipes::Dict{Int,Pipeline},
#     P::Vector{Int},
#     compressors::Dict{Int,Compressor},
#     C::Vector{Int},
#     T::Vector{Int})
#     """
#     Ensures balance between gas injections and extractions at each node i at time step t.
#     Accounts for flows through compressors if there are any.
#     """

#     q_s, q_d_cur, m⁺_in, m⁺_out, m⁻_in, m⁻_out = 
#         model[:q_s], model[:q_d_cur], model[:m⁺_in], model[:m⁺_out], model[:m⁻_in], model[:m⁻_out]

#     nodal_gas_balance = [
#         sum(q_s[s,t] for s in S if sources[s].node == i; init = 0.0) -
#         sum(demands[d].Q_d[t] for d in D if demands[d].node == i; init = 0.0) +
#         sum(q_d_cur[d,t] for d in D if demands[d].node == i; init = 0.0) +
#         sum(m_out[p,t] for p in P if pipes[p].stop == i; init = 0.0) -
#         sum(m_in[p,t] for p in P if pipes[p].start == i; init = 0.0)
#             for i in I, t in T]

#     if isempty(compressors) == true
#         @expression(model, expr_nodal_gas_balance[i in I, t in T], nodal_gas_balance[i,t])
#     else
#         m_c = model[:m_c]
#         @expression(model, expr_nodal_gas_balance[i in I, t in T], 
#             nodal_gas_balance[i,t] +
#             sum(m_c[c,t] for c in C if compressors[c].stop == i) -
#             sum(m_c[c,t] for c in C if compressors[c].start == i))
#     end
#     return nothing
# end


function constraint_nodal_gas_trade_balance(
    model::JuMP.Model,
    I::Vector{Int},
    sources::Dict{Int,GasSource},
    S::Vector{Int},
    demands::Dict{Int,GasLoad},
    D::Vector{Int},
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    compressors::Dict{Int,Compressor},
    C::Vector{Int},
    T::Vector{Int})
    """
    Ensures balance between traded gas inflows and outflows at each node i at time step t.
    """

    q_s, q_d_cur, m = 
        model[:q_s], model[:q_d_cur], model[:m]
    @expression(model, expr_nodal_gas_balance[i in I, t in T],
        sum(q_s[s,t] for s in S if sources[s].node == i) -
        sum(demands[d].Q_d[t] for d in D if demands[d].node == i) +
        sum(q_d_cur[d,t] for d in D if demands[d].node == i) +
        sum(m[p,t] for p in P if pipes[p].stop == i) -
        sum(m[p,t] for p in P if pipes[p].start == i))

    if isempty(compressors) == true
        @constraint(model, nodal_gas_balance[i in I, t in T], expr_nodal_gas_balance[i,t] == 0)
    else
        m_c = model[:m_c]
        @constraint(model, nodal_gas_balance[i in I, t in T], 
            expr_nodal_gas_balance[i,t] +
            sum(m_c[c,t] for c in C if compressors[c].stop == i) -
            sum(m_c[c,t] for c in C if compressors[c].start == i) == 0)
    end

    return nothing
end


function constraint_subpipeline_flow(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P_hat::Vector{Int},
    T::Vector{Int})
    """
    Defines the average flow m and pressure π_avg in pipeline p at timestep t.
    """
    m = model[:m]

    P_hat_sub = [[p for (p,pipe) in pipes if pipe.parent == p_hat] for p_hat in P_hat]
    no_subpipes = length.(P_hat_sub)

    @constraints model begin
        subpipeflow[p_hat in P_hat, p_tup_no in 1:no_subpipes[p_hat]-1, t in T],
            m[P_hat_sub[p_hat][p_tup_no],t] * m[P_hat_sub[p_hat][p_tup_no+1],t] >= 0
    end
    return nothing
end


function constraint_average_flow_and_pressure(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    T::Vector{Int})
    """
    Defines the average flow m and pressure π_avg in pipeline p at timestep t.
    """
    m, m_in, m_out, π, π_avg = model[:m], model[:m_in], model[:m_out], model[:π], model[:π_avg]  
    @constraints model begin
        avg_flow[p in P, t in T],
            m[p,t] == (m_in[p,t]+m_out[p,t])/2
        avg_pressure[p in P, t in T],
            π_avg[p,t] == (π[pipes[p].start,t]+π[pipes[p].stop,t])/2
    end
    return nothing
end


function constraint_conservation_mass(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    ut1::Int,
    dt::Real,
    T::Vector{Int})
    """
    Ensures the conservation of mass in pipeline p at timestep t.
    """
    m_in, m_out, π_avg = model[:m_in], model[:m_out], model[:π_avg]
    # @constraints model begin
    #     conservation_mass_initial[p in P, t=[1]],
    #         pipes[p].C^2/pipes[p].A*(m_out[p,t]-m_in[p,t])/pipes[p].L == 0
    #     conservation_mass[p in P, t in T[2:end]],
    #         ut1*(π_avg[p,t]-π_avg[p,t-1])/dt +
    #         pipes[p].C^2/pipes[p].A*(m_out[p,t]-m_in[p,t])/pipes[p].L == 0

    # end
    @constraint(model, conservation_mass[p in P, t in T],
        ut1*(π_avg[p,t]-π_avg[p,t-1])/dt +
        pipes[p].C^2/pipes[p].A*(m_out[p,t]-m_in[p,t])/pipes[p].L == 0)
    return nothing
end


function constraint_conservation_momentum(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    ut2::Int,
    dt::Real,
    T::Vector{Int})
    """
    Ensures the conservation of mass and momentum in pipeline p at timestep t.
    Utilizes auxiliary variable γ to replace the nonconvex m*|m| term.
    """
    m, π, γ = model[:m], model[:π], model[:γ]
    # @constraints model begin
    #     conservation_momentum_initial[p in P, t=[1]],
    #         pipes[p].A*(π[pipes[p].stop,t]-π[pipes[p].start,t])/pipes[p].L +
    #         pipes[p].F*pipes[p].C^2/(2*pipes[p].D*pipes[p].A)*γ[p,t] == 0
    #     conservation_momentum[p in P, t in T[2:end]],
    #         ut2*(m[p,t]-m[p,t-1])/dt +
    #         pipes[p].A*(π[pipes[p].stop,t]-π[pipes[p].start,t])/pipes[p].L +
    #         pipes[p].F*pipes[p].C^2/(2*pipes[p].D*pipes[p].A)*γ[p,t] == 0
    # end
    @constraint(model, conservation_momentum[p in P, t in T],
        ut2*(m[p,t]-m[p,t-1])/dt +
        pipes[p].A*(π[pipes[p].stop,t]-π[pipes[p].start,t])/pipes[p].L +
        pipes[p].F*pipes[p].C^2/(2*pipes[p].D*pipes[p].A)*γ[p,t] == 0)
    return nothing
end


function constraint_conservation_momentum_nonconvex_twoway(
    model::JuMP.Model,
    P::Vector{Int},
    T::Vector{Int})
    """
    Enforces the nonconvex relationship between flows and pressures in pipeline p at timestep t based on
    twoway flows. Only works with Ipopt Solver.
    """
    m, π_avg, γ = model[:m], model[:π_avg], model[:γ]
    @NLconstraint(model, conservation_momentum_nonconvex[p in P, t in T],
        (abs(m[p,t])*m[p,t]) == γ[p,t]*π_avg[p,t])
    return nothing
end


function constraint_conservation_momentum_soc(
    model::JuMP.Model,
    P::Vector{Int},
    T::Vector{Int})
    """
    Enforces the nonconvex relationship between flows and pressures in
    pipeline p at timestep t based on oneway flows.
    Based on conic relaxations. Works with every MISOCP solver.
    """
    m⁺, m⁻, π_avg, γ⁺, γ⁻ = model[:m⁺], model[:m⁻], model[:π_avg], model[:γ⁺], model[:γ⁻]
    @constraints model begin
        pos_conservation_momentum_convex[p in P, t in T],
            m⁺[p,t]*m⁺[p,t] <= π_avg[p,t]*γ⁺[p,t]
            # [γ⁺[p,t],1/2*π_avg[p,t], m⁺[p,t]] in RotatedSecondOrderCone()
        neg_conservation_momentum_convex[p in P, t in T],
            m⁻[p,t]*m⁻[p,t] <= π_avg[p,t]*γ⁻[p,t]
            # [γ⁻[p,t],1/2*π_avg[p,t], m⁻[p,t]] in RotatedSecondOrderCone()
    end
    return nothing
end


function constraint_conservation_momentum_soc_McCormick(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    T::Vector{Int})
    """
    Enforces the nonconvex relationship between flows and pressures in
    pipeline p at timestep t based on oneway flows.
    Based on conic relaxations and McCormick envelop. Works with every MISOCP solver.
    """
    m⁺, m⁻, v⁺, v⁻ = model[:m⁺], model[:m⁻], model[:v⁺], model[:v⁻]
    @constraints model begin
        pos_conservation_momentum_convex[p in P, t in T],
            m⁺[p,t]*m⁺[p,t] <= v⁺[p,t]
        neg_conservation_momentum_convex[p in P, t in T],
            m⁻[p,t]*m⁻[p,t] <= v⁻[p,t]
    end

    π_avg = model[:π_avg]  
    γ⁺, v⁺ = model[:γ⁺], model[:v⁺]
    nf_mc2 = [
        10^(floor(maximum(log10.(abs.([pipes[p].Γ_max*pipes[p].Π_avg_max, pipes[p].Γ_max])))) - 2)
            for p in P]
    nf_mc3 = [
        10^(floor(maximum(log10.(abs.([pipes[p].Γ_max*pipes[p].Π_avg_min, pipes[p].Γ_max])))) - 2)
            for p in P]

    @constraints model begin
        mccormick_1_pos[p in P, t in T],
            v⁺[p,t] >= γ⁺[p,t]*pipes[p].Π_avg_min
        mccormick_2_pos[p in P, t in T],
            1/nf_mc2[p] * v⁺[p,t] >= 1/nf_mc2[p] .*(pipes[p].Γ_max*π_avg[p,t] + γ⁺[p,t]*pipes[p].Π_avg_max - 
                pipes[p].Γ_max*pipes[p].Π_avg_max)
        mccormick_3_pos[p in P, t in T],
            1/nf_mc3[p] * v⁺[p,t] <= 1/nf_mc3[p] .* (pipes[p].Γ_max*π_avg[p,t] + γ⁺[p,t]*pipes[p].Π_avg_min - 
                pipes[p].Γ_max*pipes[p].Π_avg_min)
        mccormick_4_pos[p in P, t in T],
            v⁺[p,t] <= γ⁺[p,t]*pipes[p].Π_avg_max
    end

    γ⁻, v⁻ = model[:γ⁻], model[:v⁻]
    nf_mc2 = [
        10^(floor(maximum(log10.(abs.([pipes[p].Γ_max*pipes[p].Π_avg_max, pipes[p].Γ_min])))) - 2)
            for p in P]
    nf_mc3 = [
        10^(floor(maximum(log10.(abs.([pipes[p].Γ_max*pipes[p].Π_avg_min, pipes[p].Γ_min])))) - 2)
            for p in P]

    @constraints model begin
        mccormick_1_neg[p in P, t in T],
            v⁻[p,t] >= γ⁻[p,t]*pipes[p].Π_avg_min
        mccormick_2_neg[p in P, t in T],
            1/nf_mc2[p] * v⁻[p,t] >= 1/nf_mc2[p] .* (-pipes[p].Γ_min*π_avg[p,t] + γ⁻[p,t]*pipes[p].Π_avg_max - 
                (-pipes[p].Γ_min)*pipes[p].Π_avg_max)
        mccormick_3_neg[p in P, t in T],
            1/nf_mc3[p] * v⁻[p,t] <= 1/nf_mc3[p] .* (-pipes[p].Γ_min*π_avg[p,t] + γ⁻[p,t]*pipes[p].Π_avg_min - 
                (-pipes[p].Γ_min)*pipes[p].Π_avg_min)
        mccormick_4_neg[p in P, t in T],
            v⁻[p,t] <= γ⁻[p,t]*pipes[p].Π_avg_max
    end


    return nothing
end

# function McCormick_envelop_1(w, x, y, X_min, Y_min)
#     return w >= X_min*y + x*Y_min - X_min*Y_min
# end

# function McCormick_envelop_2(w, x, y, X_max, Y_max)
#     return w >= X_max*y + x*Y_max - X_mmax*Y_max
# end

# function McCormick_envelop_3(w, x, y, X_max, Y_min)
#     return w <= X_max*y + x*Y_min - X_max*Y_min
# end

# function McCormick_envelop_4(w, x, y, X_min, Y_max)
#     return w <= X_min*y + x*Y_max - X_min*Y_max
# end

function constraint_conservation_momentum_nonconvex_oneway(
    model::JuMP.Model,
    P::Vector{Int},
    T::Vector{Int})
    """
    Enforces the nonconvex relationship between flows and pressures in
    pipeline p at timestep t based on oneway flows.
    Only works with Gurobi Solver and NonConvex = 2.
    """
    m⁺, m⁻, π_avg, γ⁺, γ⁻ = model[:m⁺], model[:m⁻], model[:π_avg], model[:γ⁺], model[:γ⁻]
    @constraints model begin
        pos_conservation_momentum_convex[p in P, t in T],
            m⁺[p,t]*m⁺[p,t] == π_avg[p,t]*γ⁺[p,t]
        neg_conservation_momentum_convex[p in P, t in T],
            m⁻[p,t]*m⁻[p,t] == π_avg[p,t]*γ⁻[p,t]
    end
    return nothing
end


function constraint_oneway_flow_limits(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    T::Vector{Int})
    """
    Defines the twoway flow as the sum of oneway flows in pipeline p at timestep t.
    Ensures that the flow can only be non-zero when the associated binary is equal to one.
    """
    m, m⁺, m⁻, z = model[:m], model[:m⁺], model[:m⁻], model[:z]
    @constraints model begin
        def_twoway_flow[p in P, t in T],
            m[p,t] == m⁺[p,t] - m⁻[p,t]
        # oneway_flow_sos1[p in P, t in T],
        #     [m⁺[p,t],m⁻[p,t]] in SOS1()
        pos_flow_limit_upper[p in P, t in T],
            m⁺[p,t] <= pipes[p].M_max * z[p,t]
        neg_flow_limit_upper[p in P, t in T],
            m⁻[p,t] <= -pipes[p].M_min * (1-z[p,t])
    end
    return nothing
end


# function constraint_oneway_flow_limits_AS(
#     model::JuMP.Model,
#     pipes::Dict{Int,Pipeline},
#     P::Vector{Int},
#     T::Vector{Int})
#     """
#     Defines the twoway flow as the sum of oneway flows in pipeline p at timestep t.
#     Ensures that the flow can only be non-zero when the associated binary is equal to one.
#     """
#     m, m⁺, m⁻, z = model[:m], model[:m⁺], model[:m⁻], model[:z]
#     @constraints model begin
#         def_twoway_flow[p in P, t in T],
#             m[p,t] == m⁺[p,t] - m⁻[p,t]
#         # oneway_flow_sos1[p in P, t in T],
#         #     [m⁺[p,t],m⁻[p,t]] in SOS1()
#         pos_flow_limit_upper[p in P, t in T],
#             m⁺[p,t] <= pipes[p].M_max * z[p,t]
#         neg_flow_limit_upper[p in P, t in T],
#             m⁻[p,t] <= -pipes[p].M_min * (1-z[p,t])
#     end
#     return nothing
# end


function constraint_average_flow_AS(
    model::JuMP.Model,
    P::Vector{Int},
    T::Vector{Int})
    """
    Defines the average flow m and pressure π_avg in pipeline p at timestep t.
    """
    m, m_in, m_out = model[:m], model[:m_in], model[:m_out]
    # m⁺_in, m⁺_out, m⁻_in, m⁻_out = model[:m⁺_in], model[:m⁺_out], model[:m⁻_in], model[:m⁻_out]

    @constraints model begin
        avg_flow[p in P, t in T],
            m[p,t] == (m_in[p,t]+m_out[p,t])/2
        # avg_flow_pos[p in P, t in T],
        #     m[p,t] == (m⁺_in[p,t] + m⁺_out[p,t])/2
        # avg_flow_neg[p in P, t in T],
        #     m[p,t] == (m⁻_in[p,t] + m⁻_out[p,t])/2
    end
    return nothing
end


function constraint_linepack_AS(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    dt::Real,
    T::Vector{Int})
    """
    Adds linepack constraints based on explicit linepack variable LP.
    """

    # m⁺_in, m⁺_out, m⁻_in, m⁻_out, LP =
    #     model[:m⁺_in], model[:m⁺_out], model[:m⁻_in], model[:m⁻_out], model[:LP]

    m_in, m_out, ϕ⁺, LP =
        model[:m_in], model[:m_out], model[:ϕ⁺], model[:LP]

    @constraints model begin
        linepack_def[p in P, t in T],
            LP[p,t] == pipes[p].S*ϕ⁺[p,t]/2
        linepack_intertemporal_initial[p in P, t = [1]],
            m_in[p,t] - m_out[p,t] == 0
        linepack_intertemporal[p in P, t in T[2:end]],
            LP[p,t] == LP[p,t-1] +
                dt * (m_in[p,t] - m_out[p,t])
        linepack_end_condition,
            sum(LP[p,T[end]] for p in P) >= sum(LP[p,T[1]] for p in P)
    end
    return nothing
end


function constraint_McCormick_envelops_AS(
    model::JuMP.Model,
    nodes::Dict{Int, Node},
    pipes::Dict{Int, Pipeline},
    P::Vector{Int}, 
    T::Vector{Int})
    
    ϕ⁺, ϕ⁻, ψ = model[:ϕ⁺], model[:ϕ⁻], model[:ψ]

    ϕ⁺_min = [nodes[pipes[p].start].Π_min + nodes[pipes[p].stop].Π_min for p in P]
    ϕ⁺_max = [nodes[pipes[p].start].Π_max + nodes[pipes[p].stop].Π_max for p in P]
    ϕ⁻_min = [nodes[pipes[p].start].Π_min - nodes[pipes[p].stop].Π_max for p in P]
    ϕ⁻_max = [nodes[pipes[p].start].Π_max - nodes[pipes[p].stop].Π_min for p in P]

    @constraints model begin
        McCormick_1[p in P, t in T],
            ψ[p,t] >= ϕ⁺_min[p] * ϕ⁻[p,t] + ϕ⁺[p,t] * ϕ⁻_min[p] - ϕ⁺_min[p] * ϕ⁻_min[p]
        McCormick_2[p in P, t in T],
            ψ[p,t] >= ϕ⁺_max[p] * ϕ⁻[p,t] + ϕ⁺[p,t] * ϕ⁻_max[p] - ϕ⁺_max[p] * ϕ⁻_max[p]
        McCormick_3[p in P, t in T],
            ψ[p,t] <= ϕ⁺_max[p] * ϕ⁻[p,t] + ϕ⁺[p,t] * ϕ⁻_min[p] - ϕ⁺_max[p] * ϕ⁻_min[p]
        McCormick_4[p in P, t in T],
            ψ[p,t] >= ϕ⁺_min[p] * ϕ⁻[p,t] + ϕ⁺[p,t] * ϕ⁻_max[p] - ϕ⁺_min[p] * ϕ⁻_max[p]
    end
    return nothing
end


function constraint_weymouth_AS(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int}, 
    T::Vector{Int}
)
    """
    Reformulated Weymouth equation.
    """
    m, ψ, z = model[:m], model[:ψ], model[:z]

    BigM_max = [sqrt(2)*pipes[p].M_max for p in P]
    BigM_min = [sqrt(2)*pipes[p].M_min for p in P]

    @constraints model begin
        weymouth_soc_pos[p in P, t in T],
            m[p,t]^2 <= pipes[p].K^2 * ψ[p,t] + BigM_min[p]^2 * (1-z[p,t])
        weymouth_soc_neg[p in P, t in T],
            m[p,t]^2 <= -pipes[p].K^2 * ψ[p,t] + BigM_max[p]^2 * z[p,t]
    end
    return nothing
end



function constraint_oneway_conservation_momentum(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    T::Vector{Int})
    """
    Defines the twoway flow as the sum of oneway flows in pipeline p at timestep t.
    Ensures that the flow can only be non-zero when the associated binary is equal to one.
    """
    γ, γ⁺, γ⁻, z = model[:γ], model[:γ⁺], model[:γ⁻], model[:z]
    @constraints model begin
        def_twoway_momentum_auxiliary[p in P, t in T],
            γ[p,t] == γ⁺[p,t] - γ⁻[p,t]       
        # momentum_auxilary_sos1[p in P, t in T],
        #     [γ⁺[p,t],γ⁻[p,t]] in SOS1()
        pos_momentum_auxiliary_limit_upper[p in P, t in T],
            γ⁺[p,t] <= pipes[p].Γ_max * z[p,t]
        negative_momentum_auxiliary_limit_upper[p in P, t in T],
            γ⁻[p,t] <= -pipes[p].Γ_min * (1-z[p,t])
    end
    return nothing
end


# function constraint_convex_envelops_conservation_momentum(
#     model::JuMP.Model,
#     pipes::Dict{Int,Pipeline},
#     P::Vector{Int},
#     T::Vector{Int})
#     """
#     Convex envelops for the nonconvex part of the 
#     conservation of momentum constraint for pipeline p at timestep t.
#     The envelop for m*|m| is based on Mhanna, S. et al. (2021). Iterative LP-based Methods
#     for the Multiperiod Optimal Electricity and Gas Flow Problem.
#     The under- and over-estimators associated with the tangent of the curve, 1 and 3, respectively,
#     are only created iff they do not limit the feasible region of m below and above M_max and M_min, respectively.
#     This is the case if the slope of the respective constraint is larger than the other under- and over-estimator,
#     2 and 4, respectively. For example, considering 3 and 4, iff -2*M_min >= (sqrt(8)-2)M_max
#     The envelop for γ*π_avg is a McCormick envelop.
#     """
#     # Convex envelop m*|m|
#     m, u = model[:m], model[:u]
#     P_pos = [p for p in P if -2*pipes[p].M_min >= (sqrt(8)-2)*pipes[p].M_max]
#     P_neg = [p for p in P if 2*pipes[p].M_max >= (2-sqrt(8))*pipes[p].M_min]
#     δ_under = [(pipes[p].M_min^2*(3-sqrt(8)) - pipes[p].M_max^2) / ((2-sqrt(8))*pipes[p].M_min - 2*pipes[p].M_max)
#         for p in P, t in T]
#     δ_over = [(pipes[p].M_min^2 - pipes[p].M_max^2*(3-sqrt(8))) / ((sqrt(8)-2)*pipes[p].M_max + 2*pipes[p].M_min)
#         for p in P, t in T]

#     @constraints model begin
#         envelop_1[p in P, t in T],
#             u[p,t] >= m[p,t]*(2-sqrt(8))*pipes[p].M_min - pipes[p].M_min^2*(3-sqrt(8))
#         envelop_2[p in P_neg, t in T],
#             u[p,t] >= 2*m[p,t]*pipes[p].M_max - pipes[p].M_max^2
#         envelop_5[p in P_neg, t in T],
#             u[p,t] >= δ_under[p,t]*abs(δ_under[p,t]) + 2*abs(δ_under[p,t])*(m[p,t]-δ_under[p,t])
#         envelop_3[p in P, t in T],
#             u[p,t] <= m[p,t]*(sqrt(8)-2)*pipes[p].M_max + pipes[p].M_max^2*(3-sqrt(8))
#         envelop_4[p in P_pos, t in T],
#             u[p,t] <= -2*m[p,t]*pipes[p].M_min + pipes[p].M_min^2
#         envelop_6[p in P_pos, t in T],
#             u[p,t] <= δ_over[p,t]*abs(δ_over[p,t]) + 2*abs(δ_over[p,t])*(m[p,t]-δ_over[p,t])
#     end
#     # Convex envelop γ*π_avg (McCormick envelop)
#     π_avg, γ, v = model[:π_avg], model[:γ], model[:v]
#     @constraints model begin
#         mccormick_1[p in P, t in T],
#             v[p,t] >= pipes[p].Γ_min*π_avg[p,t] + γ[p,t]*pipes[p].Π_avg_min - 
#                 pipes[p].Γ_min*pipes[p].Π_avg_min
#         mccormick_2[p in P, t in T],
#             v[p,t] >= pipes[p].Γ_max*π_avg[p,t] + γ[p,t]*pipes[p].Π_avg_max - 
#                 pipes[p].Γ_max*pipes[p].Π_avg_max
#         mccormick_3[p in P, t in T],
#             v[p,t] <= pipes[p].Γ_max*π_avg[p,t] + γ[p,t]*pipes[p].Π_avg_min - 
#                 pipes[p].Γ_max*pipes[p].Π_avg_min
#         mccormick_4[p in P, t in T],
#             v[p,t] <= pipes[p].Γ_min*π_avg[p,t] + γ[p,t]*pipes[p].Π_avg_max - 
#                 pipes[p].Γ_min*pipes[p].Π_avg_max
#     end

#     # constraint_McCormick_envelop(model, pipes, P, T)

#     # Equality based on nonconvex momentum conversvation term
#     @constraint(model, envelop_auxiliary_equality[p in P, t in T], u[p,t] == v[p,t])
#     return nothing
    
# end


function constraint_convex_envelops_conservation_momentum(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    nodes::Dict{Int,<:Node},
    T::Vector{Int})
    """
    Convex envelops for the nonconvex part of the 
    conservation of momentum constraint for pipeline p at timestep t.
    This is a polyhedral envelop based on the idea in Mhanna, S. et al. (2021). Iterative LP-based Methods
    for the Multiperiod Optimal Electricity and Gas Flow Problem.
    """
    m, π_avg, γ = model[:m], model[:π_avg], model[:γ] 

    Π_diff_max_pos = [(nodes[pipes[p].start].Π_max+nodes[pipes[p].stop].Π_min)/2 for p in P]
    Π_diff_max_neg = [(nodes[pipes[p].stop].Π_max+nodes[pipes[p].start].Π_min)/2 for p in P]
    M_intersection_o = [
        (-pipes[p].M_max^2*(3-sqrt(8))+pipes[p].M_min^2)/(pipes[p].M_max*(sqrt(8)-2)+2*pipes[p].M_min) for p in P]
    M_intersection_u = [
        (-pipes[p].M_min^2*(sqrt(8)-3)-pipes[p].M_max^2)/(pipes[p].M_min*(2-sqrt(8))-2*pipes[p].M_max) for p in P]

    P_pos = [p for p in P if -2*pipes[p].M_min >= (sqrt(8)-2)*pipes[p].M_max]
    P_neg = [p for p in P if 2*pipes[p].M_max >= (2-sqrt(8))*pipes[p].M_min]

    @constraints model begin
        envelop_1[p in P, t in T],
            γ[p,t] >= m[p,t]*(2-sqrt(8))*pipes[p].M_min/Π_diff_max_pos[p] +
                (pipes[p].M_min/Π_diff_max_pos[p])^2*(sqrt(8)-3)*π_avg[p,t]
        envelop_2[p in P_neg, t in T],
            γ[p,t] >= 2*m[p,t]*pipes[p].M_max/Π_diff_max_pos[p] -
                (pipes[p].M_max/Π_diff_max_pos[p])^2*π_avg[p,t]
        envelop_3[p in P_neg, t in T],
            γ[p,t] >= 2*m[p,t]*M_intersection_u[p]/Π_diff_max_pos[p] -
                (M_intersection_u[p]/Π_diff_max_pos[p])^2*π_avg[p,t]
        envelop_4[p in P, t in T],
            γ[p,t] <= m[p,t]*(sqrt(8)-2)*pipes[p].M_max/Π_diff_max_neg[p] +
                (pipes[p].M_max/Π_diff_max_neg[p])^2*(3-sqrt(8))*π_avg[p,t]
        envelop_5[p in P_pos, t in T],
            γ[p,t] <= -2*m[p,t]*pipes[p].M_min/Π_diff_max_neg[p] +
                (pipes[p].M_min/Π_diff_max_neg[p])^2*π_avg[p,t]
        envelop_6[p in P_pos, t in T],
            γ[p,t] <= -2*m[p,t]*M_intersection_o[p]/Π_diff_max_neg[p] +
                (M_intersection_o[p]/Π_diff_max_neg[p])^2*π_avg[p,t]
    end

    return nothing
    
end


function constraint_linear_relaxation_conservation_momentum(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    nodes::Dict{Int,<:Node},
    T::Vector{Int})
    """
    Mixed-integer linear relaxation of the nonconvex part of the conservation
    of momentum equation.
    """
    m⁺, m⁻, π_avg, γ⁺, γ⁻ = model[:m⁺], model[:m⁻], model[:π_avg], model[:γ⁺], model[:γ⁻] 

    Π_diff_max_pos = [(nodes[pipes[p].start].Π_max+nodes[pipes[p].stop].Π_min)/2 for p in P]
    Π_diff_max_neg = [(nodes[pipes[p].stop].Π_max+nodes[pipes[p].start].Π_min)/2 for p in P]
    M_intersection_o_1 = [
        (-pipes[p].M_max^2*(3-sqrt(8))+pipes[p].M_min^2)/(pipes[p].M_max*(sqrt(8)-2)+2*pipes[p].M_min) for p in P]
    M_intersection_u_1 = [
        (-pipes[p].M_min^2*(sqrt(8)-3)-pipes[p].M_max^2)/(pipes[p].M_min*(2-sqrt(8))-2*pipes[p].M_max) for p in P]

    M_intersection_o_2 = [-pipes[p].M_max*(sqrt(8)-3)/(2-sqrt(8)) for p in P]
    M_intersection_u_2 = [-pipes[p].M_min*(sqrt(8)-3)/(2-sqrt(8)) for p in P]

    @constraints model begin
        envelop_1[p in P, t in T],
            γ⁺[p,t] >= m⁺[p,t]*(2-sqrt(8))*pipes[p].M_min/Π_diff_max_pos[p] +
                (pipes[p].M_min/Π_diff_max_pos[p])^2*(sqrt(8)-3)*π_avg[p,t]
        envelop_2[p in P, t in T],
            γ⁺[p,t] >= 2*m⁺[p,t]*pipes[p].M_max/Π_diff_max_pos[p] -
                (pipes[p].M_max/Π_diff_max_pos[p])^2*π_avg[p,t]
        envelop_3[p in P, t in T],
            γ⁺[p,t] >= 2*m⁺[p,t]*M_intersection_u_1[p]/Π_diff_max_pos[p] -
                (M_intersection_u_1[p]/Π_diff_max_pos[p])^2*π_avg[p,t]
        envelop_7[p in P, t in T],
            γ⁺[p,t] >= 2*m⁺[p,t]*M_intersection_u_2[p]/Π_diff_max_pos[p] -
                (M_intersection_u_2[p]/Π_diff_max_pos[p])^2*π_avg[p,t]
        envelop_4[p in P, t in T],
            γ⁻[p,t] >= m⁻[p,t]*(sqrt(8)-2)*pipes[p].M_max/Π_diff_max_neg[p] -
                (pipes[p].M_max/Π_diff_max_neg[p])^2*(3-sqrt(8))*π_avg[p,t]
        envelop_5[p in P, t in T],
            γ⁻[p,t] >= -2*m⁻[p,t]*pipes[p].M_min/Π_diff_max_neg[p] -
                (pipes[p].M_min/Π_diff_max_neg[p])^2*π_avg[p,t]
        envelop_6[p in P, t in T],
            γ⁻[p,t] >= -2*m⁻[p,t]*M_intersection_o_1[p]/Π_diff_max_neg[p] -
                (M_intersection_o_1[p]/Π_diff_max_neg[p])^2*π_avg[p,t]
        envelop_8[p in P, t in T],
            γ⁻[p,t] >= -2*m⁻[p,t]*M_intersection_o_2[p]/Π_diff_max_neg[p] -
                (M_intersection_o_2[p]/Π_diff_max_neg[p])^2*π_avg[p,t]
    end

    return nothing
    
end


function constraint_McCormick_envelop(
    w, x, y, X_min, X_max, Y_min, Y_max, P, T
    )
    """
    w = x*y

    w >= X_min*y + x*Y_min - X_min*Y_min
    w >= X_max*y + x*Y_max - X_mmax*Y_max
    w <= X_max*y + x*Y_min - X_max*Y_min
    w <= X_min*y + x*Y_max - X_min*Y_max
    """

    # π_avg, γ, v = model[:π_avg], model[:γ], model[:v]
    @constraints model begin
        mccormick_1[p in P, t in T],
            w[p,t] >= X_min[p]*y[p,t] + x[p,t]*Y_min[p] - X_min[p]*Y_min[p]
        mccormick_2[p in P, t in T],
            w[p,t] >= X_max[p]*y[p,t] + x[p,t]*Y_max[p] - X_max[p]*Y_max[p]
        mccormick_3[p in P, t in T],
            w[p,t] <= X_max[p]*y[p,t] + x[p,t]*Y_min[p] - X_max[p]*Y_min[p]
        mccormick_4[p in P, t in T],
            w[p,t] <= X_min[p]*y[p,t] + x[p,t]*Y_max[p] - X_min[p]*Y_max[p]
    end
    return nothing
end
    
    
function constraint_PWL_ldcc(
    model::JuMP.Model,
    no_bits::Int,
    no_vertices::Int,
    simplices_mat::Matrix{Dict},
    simplices_ids::Vector{String},
    P::Vector{Int},
    T::Vector{Int})
    """
    Piecewise-linearization based on the logarithmic disaggregated convex combination method from
    Geißler, B. et al. (2012). Using piecewise linear functions for solving MINLPs.
    """

    m, π_avg, γ, λ, u = model[:m], model[:π_avg], model[:γ], model[:λ], model[:u]
    # Data tuple: (m, π_avg, γ)
    @constraints model begin
        def_m[p in P, t in T],
            m[p,t] ==  sum(
                sum(λ[p,t,s_id,j]*simplices_mat[p,t][s_id][j][1] for j in 1:no_vertices)
                for s_id in simplices_ids)
        def_π_avg[p in P, t in T],
            π_avg[p,t] ==  sum(
                sum(λ[p,t,s_id,j]*simplices_mat[p,t][s_id][j][2] for j in 1:no_vertices)
                for s_id in simplices_ids)
        def_γ[p in P, t in T],
            γ[p,t] ==  sum(
                sum(λ[p,t,s_id,j]*simplices_mat[p,t][s_id][j][3] for j in 1:no_vertices)
                for s_id in simplices_ids)
    end

    @constraints model begin
        λ_integrality[p in P, t in T],
            sum(sum(λ[p,t,s_id,j] for j in 1:no_vertices)
                for s_id in simplices_ids) == 1
        λ_active[p in P, t in T, l in 1:no_bits],
            sum(sum(parse(Int, s_id[l])*λ[p,t,s_id,j] for j in 1:no_vertices) 
                for s_id in simplices_ids) <= u[p,t,l]
        λ_inactive[p in P, t in T, l in 1:no_bits],
            sum(sum((1-parse(Int, s_id[l]))*λ[p,t,s_id,j] for j in 1:no_vertices) 
                for s_id in simplices_ids) <= 1 - u[p,t,l]
    end
    return nothing
end


function constraint_linear_overestimator(
    model::JuMP.Model,
    pipes::Dict{Int,Pipeline},
    P::Vector{Int},
    nodes::Dict{Int,<:Node},
    T::Vector{Int}
    )
    """
    Add a linear overestimator to the nonconvex part of the conservation of momentum equation.
    """

    m⁺, m⁻, γ⁺, γ⁻ = model[:m⁺], model[:m⁻], model[:γ⁺], model[:γ⁻]

    Π_diff_max_pos = [(nodes[pipes[p].start].Π_max+nodes[pipes[p].stop].Π_min)/2 for p in P]
    Π_diff_max_neg = [(nodes[pipes[p].stop].Π_max+nodes[pipes[p].start].Π_min)/2 for p in P]

    @constraints model begin
        linear_overestimator_pos[p in P, t in T],
            γ⁺[p,t] <= m⁺[p,t] * (pipes[p].M_max/Π_diff_max_pos[p])
        linear_overestimator_neg[p in P, t in T],
            γ⁻[p,t] <= m⁻[p,t] * (-pipes[p].M_min/Π_diff_max_neg[p])
    end

    return nothing
end


#################### Power system constraints ####################
function constraint_power_balance_angles(
    model::JuMP.Model,
    N::Vector{Int},
    generators::Dict{Int,DispatchableGenerator},
    G::Vector{Int}, 
    windgenerators::Dict{Int,WindGenerator},
    J::Vector{Int}, 
    demands_EL::Dict{Int,PowerDemand},
    D_el::Vector{Int}, 
    lines::Dict{Int,Line},
    L::Vector{Int},
    T::Vector{Int})
    """
    Ensures balance between power injections and withdrawls at each node n at time step t.
    This formulation is based on the bus voltage angles.
    """
    p, w, θ, p_d_cur = model[:p], model[:w], model[:θ], model[:p_d_cur]
    @constraint(model, bus_power_balance[n in N, t in T],
        sum(p[g,t] for g in G if generators[g].bus == n) + 
        sum(w[j,t] for j in J if windgenerators[j].bus == n) - 
        sum(demands_EL[d].D_el[t]-p_d_cur[d,t] for d in D_el if demands_EL[d].bus == n) - 
        sum(1/lines[l].X*(θ[lines[l].start,t]-θ[lines[l].stop,t])
            for l in L if lines[l].start == n) +
        sum(1/lines[l].X*(θ[lines[l].start,t]-θ[lines[l].stop,t])
            for l in L if lines[l].stop == n) == 0)
    return nothing
end


function constraint_power_flow_limits_angles(
    model::JuMP.Model, lines::Dict{Int,Line}, L::Vector{Int}, T::Vector{Int})
    """
    Ensures that the power flow is within line limits for each line l at timestep t.
    """
    θ = model[:θ]
    @constraints model begin
        power_flow_limit_up[l in L, t in T],
            1/lines[l].X*(θ[lines[l].start,t]-θ[lines[l].stop,t]) <= lines[l].F_max
        power_flow_limit_down[l in L, t in T],
            1/lines[l].X*(θ[lines[l].start,t]-θ[lines[l].stop,t]) >= -lines[l].F_max
    end
    return nothing
end


function constraint_voltage_angle_slack(
    model::JuMP.Model, buses::Dict{Int,Bus}, T::Vector{Int})
    """
    Fixes the voltage angle at the reference bus n_ref to zero at each time step t.
    """
    n_ref = findall(x->x.slack==1, buses)[1]
    θ = model[:θ]
    @constraint(model, voltage_angle_ref_bus[n=[n_ref], t in T], θ[n,t] == 0)
    return nothing
end

# function constraint_power_balance_ptdf(
#     model::JuMP.Model,
#     G::Vector{Int}, 
#     J::Vector{Int}, 
#     demands_EL::Vector{PowerDemand},
#     D_el::Vector{Int}, 
#     T::Vector{Int})
#     """
#     Ensures balance between power injections and withdrawls at each node n at time step t.
#     This formulation is based on the Power Transfer Distribution Factors (PTDFs).
#     """
#     @constraint(model, bus_power_balance[t in T],
#         sum(p[g,t] for g in G) + 
#         sum(w[j,t] for j in J) - 
#         sum(demands_EL[d].D_el[t] for d in D_el) == 0)
#     return nothing
# end


# function constraint_power_flow_limits_ptdf(
#     model::JuMP.Model,
#     generators::Vector{DispatchableGenerator},
#     G::Vector{Int}, 
#     windgenerators::Vector{WindGenerator},
#     J::Vector{Int}, 
#     demands_EL::Vector{PowerDemand},
#     D_el::Vector{Int}, 
#     lines::Vector{Lines},
#     L::Vector{Int},
#     T::Vector{Int})
#     """
#     Ensures that the power flow is within line limits for each line l at timestep t.
#     """
#     @constraints model begin
#         power_flow_limit_up[l in L, t in T],
#             lines[l].ptdf[[generators[g].bus for g in G]]'*([p[g, t] for g in G]) +
#             lines[l].ptdf[[windgenerators[j].bus for j in J]]'*([w[j, t] for j in J]) -
#             lines[l].ptdf[[demands_EL[d].bus for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) <= lines[l].F_max
#         power_flow_limit_down[l in L, t in T],
#             lines[l].ptdf[[generators[g].bus for g in G]]'*([p[g, t] for g in G]) +
#             lines[l].ptdf[[windgenerators[j].bus for j in J]]'*([w[j, t] for j in J]) -
#             lines[l].ptdf[[demands_EL[d].bus for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) >= -lines[l].F_max
#     end
#     return nothing
# end


# function constraint_ramping(
#     model::JuMP.Model,
#     generators::Vector{DispatchableGenerator},
#     G::Vector{Int}, 
#     dt::Real,
#     T::Vector{Int})
#     @constraints model begin
#         ramping_up[t in T[2:end], g in G], 
#             p[g,t]-p[g,t-1] <= generators[g].P_up/3600*dt
#         ramping_down[t in T[2:end], g in G],
#             p[g,t]-p[g,t-1] >= -generators[g].P_down/3600*dt
#     end
#     return nothing
# end