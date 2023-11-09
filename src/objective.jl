const VARIABLE_CONTAINER_TYPE = JuMP.Containers.DenseAxisArray{
    VariableRef,
    2,
    Tuple{Vector{Int64}, Vector{Int64}},
    Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}},
    JuMP.Containers._AxisLookup{Dict{Int64, Int64}}}}


function cost_gas_linear(
    dt::Number,
    q_s::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    sources::Dict{Int, GasSource},
    q_d_cur::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    C_cur_load_gas::Number,
    S::Vector{Int},
    D::Vector{Int},
    T::Vector{Int})
    return dt/3600 *(
        sum(q_s[s,t]*sources[s].C_1 for s in S, t in T) +
        sum(q_d_cur[d,t]*C_cur_load_gas for d in D, t in T))
end


function cost_gas_quadratic(
    dt::Number,
    q_s::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    sources::Dict{Int, GasSource},
    S::Vector{Int},
    T::Vector{Int})
    return dt/3600 *(sum(q_s[s,t]^2*sources[s].C_2 for s in S, t in T))
end


function cost_gas_compression(
    π::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    T::Vector{Int})
    return sum(compressors[c].C_compr*(π[compressors[c].stop,t]-π[compressors[c].start,t])
        for c in C, t in T; init=0)
end


function cost_gas_compression_AS(
    ϕ⁻::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    T::Vector{Int})
    return sum(compressors[c].C_compr*(-ϕ⁻[c,t]) for c in C, t in T; init=0)
end


function total_cost_gas(
    objective_type::String,
    dt::Number,
    q_s::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    sources::Dict{Int, GasSource},
    q_d_cur::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    C_cur_load_gas::Number,
    π::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    S::Vector{Int},
    D::Vector{Int},
    T::Vector{Int})

    cost_gas = cost_gas_linear(
        dt, q_s, sources, q_d_cur, C_cur_load_gas, S, D, T)
    cost_gas += objective_type == "quadratic" ? 
        cost_gas_quadratic(dt, q_s, sources, S, T) : 0
    # cost_gas += isempty(compressors)==true ? 
    #     cost_gas_compression(π, compressors, C, T) : 0

    return cost_gas
end


function total_cost_gas_AS(
    objective_type::String,
    dt::Number,
    q_s::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    sources::Dict{Int, GasSource},
    q_d_cur::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    C_cur_load_gas::Number,
    ϕ⁻::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    compressors::Dict{Int, Compressor},
    C::Vector{Int},
    S::Vector{Int},
    D::Vector{Int},
    T::Vector{Int})

    cost_gas = cost_gas_linear(
        dt, q_s, sources, q_d_cur, C_cur_load_gas, S, D, T)
    cost_gas += objective_type == "quadratic" ? 
        cost_gas_quadratic(dt, q_s, sources, S, T) : 0
    # cost_gas += isempty(compressors)==true ? 
    #     cost_gas_compression_AS(ϕ⁻, compressors, C, T) : 0

    return cost_gas
end


function expression_cost_gas_system(
    model_type::String,
    objective_type::String,
    model::JuMP.Model,
    sources::Dict{Int,GasSource},
    S::Vector{Int},
    compressors::Dict{Int,Compressor},
    C::Vector{Int},
    P::Vector{Int},
    D::Vector{Int},
    C_cur_load_gas::Real,
    dt::Real,
    T::Vector{Int})
    """
    Total gas system cost.
    """
    @assert objective_type in ["linear","quadratic"]

    q_s, q_d_cur = model[:q_s], model[:q_d_cur]
    if ~(model_type == "MISOCP_weymouth")
        π = model[:π]
        cost_gas = total_cost_gas(
            objective_type, dt, q_s, sources, q_d_cur, C_cur_load_gas, π, compressors, C, S, D, T)
    elseif model_type == "MISOCP_weymouth"
        ϕ⁻ = model[:ϕ⁻]
        cost_gas = total_cost_gas_AS(
            objective_type, dt, q_s, sources, q_d_cur, C_cur_load_gas, ϕ⁻, compressors, C, S, D, T)
    end

    @expression(model, cost_gas_system, cost_gas)
    return nothing
end


function cost_power_linear(
    dt::Number,
    p::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    generators::Dict{Int,DispatchableGenerator},
    nonGFPPs::Vector{Int},
    p_d_cur::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    D_el::Vector{Int},
    C_cur_load_power::Number,
    T::Vector{Int})
    return dt/3600 *(
        sum(generators[g].C_l*p[g,t] for g in nonGFPPs, t in T)+
        sum(p_d_cur[d,t]*C_cur_load_power for d in D_el, t in T))
end


function cost_power_quadratic(
    dt::Number,
    p::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    generators::Dict{Int,DispatchableGenerator},
    nonGFPPs::Vector{Int},
    T::Vector{Int})
    return dt/3600 *(
        sum(generators[g].C_q*p[g,t]^2 for g in nonGFPPs, t in T))
end


function total_cost_power(
    objective_type::String,
    dt::Number,
    p::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    generators::Dict{Int,DispatchableGenerator},
    nonGFPPs::Vector{Int},
    p_d_cur::Union{VARIABLE_CONTAINER_TYPE, Matrix{<:Number}},
    D_el::Vector{Int},
    C_cur_load_power::Number,
    T::Vector{Int})

    cost_power = cost_power_linear(
        dt, p, generators, nonGFPPs, p_d_cur, D_el, C_cur_load_power, T)
    cost_power += objective_type == "quadratic" ? 
        cost_power_quadratic(dt, p, generators, nonGFPPs, T) : 0

    return cost_power
end



function expression_cost_power_system(
    objective_type::String,
    model::JuMP.Model,
    generators::Dict{Int,DispatchableGenerator},
    nonGFPPs::Vector{Int},
    D_el::Vector{Int},
    C_cur_load_power::Number,
    dt::Real,
    T::Vector{Int})
    """
    Total gas system cost.
    """
    p, p_d_cur = model[:p], model[:p_d_cur]
    @expression(model, cost_power_system, total_cost_power(
        objective_type, dt, p, generators, nonGFPPs, p_d_cur, D_el, C_cur_load_power, T))
    return nothing
end


function objective_cost_gas_system(
    model_type::String,
    objective_type::String,
    model::JuMP.Model,
    sources::Dict{Int,GasSource},
    S::Vector{Int},
    compressors::Dict{Int,Compressor},
    C::Vector{Int},
    P::Vector{Int},
    D::Vector{Int},
    C_cur_load_gas::Real,
    dt::Real,
    T::Vector{Int})
    """
    Minimization of total gas system cost.
    """

    expression_cost_gas_system(
        model_type, objective_type, model, sources, S, compressors, C, P, D, C_cur_load_gas, dt, T)
    @objective(model, Min, model[:cost_gas_system])
    return nothing
end


function objective_cost_integrated_energy_system(
    model_type::String,
    objective_type::String,
    model::JuMP.Model,
    sources::Dict{Int,GasSource},
    S::Vector{Int},
    compressors::Dict{Int,Compressor},
    C::Vector{Int},
    P::Vector{Int},
    generators::Dict{Int,DispatchableGenerator},
    nonGFPPs::Vector{Int},
    D::Vector{Int},
    C_cur_load_gas::Real,
    D_el::Vector{Int},
    C_cur_load_power::Number,
    dt::Real,
    T::Vector{Int})
    """
    Minimization of total power and gas system cost.
    """

    expression_cost_gas_system(
        model_type, objective_type, model, sources, S, compressors, C, P, D, C_cur_load_gas, dt, T)
    expression_cost_power_system(objective_type, model, generators, nonGFPPs, D_el, C_cur_load_power, dt, T)
    @objective(model, Min, model[:cost_gas_system]+model[:cost_power_system])
    return nothing
end


# function objective_cost_minimization(GS::GasSystem)
#     objective_cost_gas_system(model, sources, S, compressors, C, D, C_cur_load_gas, dt, T)
# end


# function objective_cost_minimization(IES::IntegratedEnergySystem)
#     objective_cost_integrated_energy_system(
#         model, sources, S, compressors, C, generators, nonGFPPs, D, C_cur_load_gas, dt, T)
# end