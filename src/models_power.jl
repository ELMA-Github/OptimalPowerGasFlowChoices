function PS_DC_model(
    model::JuMP.Model,
    buses::Dict{Int,Bus},
    N::Vector{},
    generators::Dict{Int,DispatchableGenerator},
    G::Vector{},
    GFPPs::Vector{},
    windgenerators::Dict{Int,WindGenerator},
    J::Vector{}, 
    demands_EL::Dict{Int,PowerDemand},
    D_el::Vector{}, 
    lines::Dict{Int,Line},
    L::Vector{},
    T::Vector{Int})

    variable_power_production_dispatchables(model, generators, G, T)
    variable_power_production_renewables(model, windgenerators, J, T)
    variable_bus_voltage_angles(model, N, T)
    variable_gas_consumption_gfpp(model, GFPPs, T)
    variable_power_curtailment(model, demands_EL, D_el, T)

    constraint_power_balance_angles(
        model, N, generators, G, windgenerators, J, demands_EL, D_el, lines, L, T)
    constraint_power_flow_limits_angles(model, lines, L, T)
    constraint_voltage_angle_slack(model, buses, T)

    return nothing
end
    

# function PS_DC_model(
#     model::JuMP.Model,
#     generators::Vector{DispatchableGenerator},
#     G::Vector{}, 
#     windgenerators::Vector{WindGenerator},
#     J::Vector{}, 
#     demands_EL::Vector{PowerDemand},
#     D_el::Vector{}, 
#     lines::Vector{Lines},
#     L::Vector{},
#     T::Vector{Int})
# 	"""
#     p: Power production of dispatchable unit g at timestep t [MW].
#     w: Power production of renewable unit j at timestep t [MW].
#     q_g: Gas consumption of gas-fired power plant g at timestep t [kg/s].
#     θ: Voltage phase angle of bus n at timestep t [rad].

#     bus_power_balance: Ensures balance between power injections and withdrawls at each
#         node n at time step t. This formulation is based on the bus voltage angles.
#     power_flow_limit_up/power_flow_limit_down: Ensures that the power flow is within
#         line limits for each line l at timestep t.
#     """

#     @variable model begin
#         0 <= p[g in G, t in T] <= generators[g].P_max
#         0 <= w[j in J, t in T] <= windgenerators[j].P_w[t]
#         0 <= q_g[g in GFPPs, t in T]
#         θ[n in N, t in T]
#     end

#     @constraints model begin
#         bus_power_balance[t in T, n in N],
#             sum(p[g,t] for g in G if generators[g].bus == n) + 
#             sum(w[j,t] for j in J if windgenerators[j].bus == n) - 
#             sum(demands_EL[d].D_el[t] for d in D_el if demands_EL[d].bus == n) - 
#             sum(1/lines[l].X*(θ[lines[l].start,t]-θ[lines[l].stop,t])
#                 for l in L if lines[l].start == n) +
#             sum(1/lines[l].X*(θ[lines[l].start,t]-θ[lines[l].stop,t])
#                 for l in L if lines[l].stop == n) == 0
#         power_flow_limit_up[t in T, l in L],
#             1/lines[l].X*(θ[lines[l].start,t]-θ[lines[l].stop,t]) <= lines[l].F_max
#         power_flow_limit_down[l in L, t in T],
#             1/lines[l].X*(θ[lines[l].start,t]-θ[lines[l].stop,t]) >= -lines[l].F_max
#     end
#     return nothing
# end