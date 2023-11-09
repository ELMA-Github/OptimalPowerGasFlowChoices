include("models_gas.jl")
include("models_power.jl")

function run_model!(
    ES::EnergySystem,
    config_dict_algorithm::Dict{Symbol,<:Any};
    config_dict_solver::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

    model_type = ES.model_type

    if ~isempty(ES.warmstart)
        warmstart!(ES, ES.warmstart[:method])
    end

    if ~(model_type in ["SCP", "MISOCP_sbnb"])
        model = Model()
        build_model(model, ES, config_dict_algorithm; model_type)
        if ~isempty(ES.warmstart)
            set_warmstart!(model, model_type, ES)
        end
        ES.model = model
        set_solver(model, model_type, config_dict_solver)
        solve_model!(model)
        write_solution!(ES, model, model_type)
    elseif model_type == "SCP"
        SCP(ES, config_dict_solver; config_dict_algorithm...)
    elseif model_type == "MISOCP_sbnb"
        SBNB(ES, config_dict_solver, config_dict_algorithm)
    end
    return nothing
end


function build_model(
    model::JuMP.Model,
    GS::GasSystem,
    config_dict_algorithm::Dict{Symbol,<:Any}=Dict{Symbol,Any}();
    model_type=nothing)

    nodes, sources, demands, pipes, compressors,
        I, S, D, P, P_hat, T, C, boundary_conditions, dt, C_cur_load_gas = 
        get_gas_system_attributes(GS)

    dt_pu = dt/GS.bases_dict[:T_base]
    
    build_gas_model(
        model_type, model, nodes, I, sources, S, demands, D, pipes, P, P_hat,
        compressors, C, boundary_conditions, dt_pu, T;
        config_dict_algorithm)

    constraint_nodal_gas_balance_GS(model, I, T)
    # constraint_linepack_end_condition_GS(model, demands, D, S, T)

    objective_cost_gas_system(
        model_type, ES.objective_type, model, sources, S, compressors, C, P, D, C_cur_load_gas, dt, T)

    if GS.print_model == 1 print(model) end

    return nothing
end


function build_model(
    model::JuMP.Model,
    IES::IntegratedEnergySystem,
    config_dict_algorithm::Dict{Symbol,<:Any}=Dict{Symbol,Any}();
    model_type=nothing)

    nodes, sources, demands, pipes, compressors,
        I, S, D, P, P_hat, T, C, boundary_conditions, dt, C_cur_load_gas =
        get_gas_system_attributes(IES)

    generators, windgenerators, demands_EL, lines,
        N, G, J, G, D_el, L, buses, C_cur_load_power = get_power_system_attributes(IES)
    GFPPs = [g for g in G if isa(generators[g], GasFiredGenerator)]
    nonGFPPs = [g for g in G if isa(generators[g], NonGasFiredGenerator)]

    dt_pu = dt/IES.bases_dict[:T_base]

    build_gas_model(
        model_type, model, nodes, I, sources, S, demands, D, pipes, P, P_hat,
        compressors, C, boundary_conditions, dt_pu, T;
        config_dict_algorithm)

    PS_DC_model(
        model, buses, N, generators, G, GFPPs, windgenerators, J, demands_EL, D_el, lines, L, T)

    constraint_nodal_gas_balance_IES(model, I, generators, G, T)


    objective_cost_integrated_energy_system(
        model_type, IES.objective_type, model, sources, S, compressors, C, P,
        generators, nonGFPPs, D, C_cur_load_gas, D_el, C_cur_load_power, dt, T)
    
    if IES.print_model == 1 print(model) end

    return nothing
end


function solve_model!(m::JuMP.Model)
    @time optimize!(m)
    status = termination_status(m)
    println(status)
    println(raw_status(m))
end


function create_default_solver(model_type::String)
    config_dict_solver = Dict{Symbol, Any}()
    if model_type in ["NCNLP", "boundary_conditions"]
        config_dict_solver[:solver_name] = "Ipopt"
        config_dict_solver[:max_cpu_time] = 7200.0
    elseif model_type in ["SCP", "MINLP", "MISOCP", "MISOCP_McCormick", "MISOCP_weymouth", "PWL", "trade", "CELP", "MILP"]
        config_dict_solver[:solver_name] = "Gurobi"
        config_dict_solver[:TimeLimit] = 1000000.0
    elseif model_type in ["MISOCP_sbnb"]
        config_dict_solver[:solver_name] = "CPLEX"
        config_dict_solver[:CPX_PARAM_TILIM] = 1000.0
    end
    return config_dict_solver
end


function set_solver(
    model::JuMP.Model,
    model_type::String,
    attribute_dict::Dict{Symbol,Any}=Dict{Symbol,Any}()
    )
    """
    Sets solver with attributes for the optimization model.
    """

    if isempty(attribute_dict)
        attribute_dict = create_default_solver(model_type)
    end

    if attribute_dict[:solver_name] == "Ipopt"
        set_optimizer(model, Ipopt.Optimizer)
        set_optimizer_attribute(model, "max_iter", 10000)
        
    elseif attribute_dict[:solver_name] == "Gurobi"
        set_optimizer(model, Gurobi.Optimizer)
        # set_optimizer_attribute(model, "Threads", 12)
        if model_type in ["MINLP"]
            @info "Setting Gurobi attribute NonConvex = 2 to solve a nonconvex model."
            set_optimizer_attribute(model, "NonConvex", 2)
        end
        if model_type in ["SCP"]
            set_optimizer_attribute(model, "OutputFlag", 1)
            set_optimizer_attribute(model, "NumericFocus", 2)
        end
    elseif attribute_dict[:solver_name] == "CPLEX"
        set_optimizer(model, CPLEX.Optimizer)
        if model_type in ["MISOCP_sbnb"]
            @info "Setting CPLEX attribute Threads = 1 and PREIND (presolve) = 0
                to use custom branching rules."
            set_optimizer_attribute(model, "CPXPARAM_Threads", 1) 
            set_optimizer_attribute(model, "CPX_PARAM_PREIND", 0)
        end
    else
        throw(ArgumentError("Solver must be either 'Ipopt', 'Gurobi', 'CPLEX'."))
    end

    delete!(attribute_dict, :solver_name)
    for (key, val) in attribute_dict
        set_optimizer_attribute(model, String(key), val)
    end

    return nothing
end


function warmstart!(ES::EnergySystem, warmstart_type::String="CELP")
    model = JuMP.Model()
    build_model(model, ES; model_type=warmstart_type)
    set_solver(model, warmstart_type)
    solve_model!(model)
    write_solution_warmstart!(ES, model, warmstart_type)
    return nothing
end


function write_solution_warmstart!(
    ES::EnergySystem, model::JuMP.Model, warmstart_type::String="CELP")
    """
    Reads specific solutions from a model, which are later used for warmstarting
    a more complex model. Estimates variable bounds for the more complex model based
    on those solutions.
    """
    
    values = Dict{Symbol, Matrix{Float64}}()

    values[:m] = [value(model[:m][p,t]) for p in ES.P, t in ES.T]
    values[:π] = [value(model[:π][i,t]) for i in ES.I, t in ES.T]
    values[:π_avg] = [value(model[:π_avg][p,t]) for p in ES.P, t in ES.T]
    values[:z] = get_flow_directions(values[:m])

    bounds = Dict{Symbol, Matrix{Float64}}()
    γ_recovered = calc_γ.(values[:m], values[:π_avg])

    bounds_expansion = 2
    # bounds[:m_upper] = abs.(values[:m])*bounds_expansion
    # bounds[:m_lower] = -abs.(values[:m])*bounds_expansion
    # bounds[:γ_upper] = abs.(γ_recovered)*bounds_expansion
    # bounds[:γ_lower] = -abs.(γ_recovered)*bounds_expansion
    bounds[:m] = abs.(values[:m])*bounds_expansion
    bounds[:γ] = abs.(γ_recovered)*bounds_expansion

    ES.warmstart[:values] = values
    ES.warmstart[:bounds] = bounds
    return nothing
end


function set_warmstart!(model::JuMP.Model, model_type::String, ES::EnergySystem)
    val, bound = ES.warmstart[:values], ES.warmstart[:bounds]
    # if model_type in ["NCLP", "MISOCP", "MINLP", "MISOCP_sbnb"]
    m, γ, π = model[:m], model[:γ], model[:π]
    z = model[:z]
    # m⁺, m⁻ = model[:m⁺], model[:m⁻]
    # γ⁺, γ⁻ = model[:γ⁺], model[:γ⁻]
    for t in ES.T
        for p in ES.P 
            set_start_value(z[p,t], val[:z][p,t])
            # set_start_value(m[p,t], val[:m][p,t])
            # set_start_value(π[i,t], val[:π][i,t])

            # set_upper_bound(m⁺[p,t], bound[:m][p,t])
            # set_upper_bound(m⁻[p,t], bound[:m][p,t])

            # set_upper_bound(γ⁺[p,t], bound[:γ][p,t])
            # set_upper_bound(γ⁻[p,t], bound[:γ][p,t])
        end
        for i in ES.I
            set_start_value(π[i,t], val[:π][i,t])
        end
     end
    # end
    return nothing
end


function get_boundary_conditions(config_dict, PDE_type)
    @info "Start deriving boundary conditions."
    system, case_study, dt = config_dict[:system], config_dict[:case_study], config_dict[:dt]
    path = joinpath(dirname(@__DIR__), "inputs", system, case_study)
    path = system == "integrated_power_and_gas" ? joinpath(path, "gas") : path
    # if system == "integrated_power_and_gas" joinpath(path, "gas") end
    segment_string = config_dict[:space_disc] == 0 ? "" : "_dx_"*string(config_dict[:segment])
    path_file = joinpath(path, "boundary_conditions_dt_$(dt)$(segment_string).csv")

    bc_df = DataFrame()
    try 
        bc_df = CSV.read(path_file, DataFrame, delim=",")
    catch ArgumentError
        gas_params = CSV.read(joinpath(path, "gas_params.csv"), DataFrame, delim=",")
        # Overwrite parameters to cover full time horizon
        config_dict[:timehorizon] = gas_params[1, :T_gasload_h]
        config_dict[:model_type] = "boundary_conditions"
        # config_dict[:dt] = gas_params[1, :dt_gasload_s]
        ES = intialize_energy_system(config_dict)
        intialize_data!(ES)

        attribute_dict = Dict(:solver_name => "Ipopt", :max_cpu_time => 7200.0)
        # First run with steady-state in initial timestep
        model = Model(Ipopt.Optimizer)
        build_model(model, ES, Dict(:PDE_type => PDE_type); model_type="boundary_conditions")
        set_solver(model, "NCNLP", copy(attribute_dict))
        optimize!(model)

        ES.boundary_conditions = Dict{Symbol,Vector}(
            :m => value.(model[:m]).data[:,end],
            :π_avg => value.(model[:π_avg]).data[:,end])

        # Second model run with boundary conditions
        model = Model(Ipopt.Optimizer)
        build_model(model, ES, Dict(:PDE_type => PDE_type); model_type="NCNLP")
        set_solver(model, "NCNLP", copy(attribute_dict))
        optimize!(model)

        # Get final boundary conditions
        bc_df = DataFrame(
            :Pipe_No => ES.P,
            :M0 => ES.bases_dict[:M_base]*value.(model[:m]).data[:,end],
            :Pavg0 => ES.bases_dict[:Π_base]*value.(model[:π_avg]).data[:,end]
        )
        # CSV.write(path_file, bc_df)
    end
    bc_dict = Dict{Symbol,Vector}(:m => bc_df[!,:M0], :π_avg => bc_df[!,:Pavg0])
    @info "Boundary conditions sucessfully derived."
    return bc_dict
end

function write_solution_vars!(ES::EnergySystem, model::JuMP.Model, vars::Vector{Symbol})
    for var in vars
        Base.setfield!(ES, var, value.(model[var]).data)
    end
    return nothing
end

function write_solution_gas!(ES::EnergySystem, model::JuMP.Model, model_type::String)

    # Model statistics
    ES.obj_val = objective_value(model)
    ES.solve_time =  JuMP.solve_time(model)

    # Common variables in all models
    write_solution_vars!(ES, model, [:q_s, :q_d_cur])
    ES.m_c = isempty(ES.compressors)==true ? [] : value.(model[:m_c]).data
    ES.m = [value(model[:m][p,t]) for p in ES.P, t in ES.T]
    # check_possible_flow_bound_violation(ES)

    if ~(model_type in ["trade", "MISOCP_weymouth"])
        write_solution_vars!(ES, model, [:π, :m_in, :m_out, :γ, :π_avg])
        # ES.π_avg = [value(model[:π_avg][p,t]) for p in ES.P, t in ES.T]
        ES.LP = calc_linepack_expost(ES)
        # Model-specific variables
        if model_type in ["MINLP", "MISOCP", "MILP"]
            write_solution_vars!(ES, model, [:m⁺, :m⁻, :γ⁺, :γ⁻])
            setfield!(ES, :z, Int.(round.(value.(model[:z]).data)))
        elseif model_type in ["NCNLP", "CELP"]
            ES.z = get_flow_directions(ES.m)
        end
    elseif model_type in ["trade"]
        ES.z = get_flow_directions(ES.m)
    elseif model_type in ["MISOCP_weymouth"]
        write_solution_vars!(ES, model, [:m_in, :m_out, :z, :LP])
        ES.π = calc_pressurs_from_ϕ(
            value.(model[:ϕ⁺]).data, value.(model[:ϕ⁻]).data,
            value.(model[:ϕ⁺_c]).data, value.(model[:ϕ⁻_c]).data,
            ES.pipes, ES.P, ES.I, ES.compressors, ES.C)
    end
    ES.total_cost = total_cost_gas(
        ES.objective_type, ES.dt, ES.q_s, ES.sources, ES.q_d_cur, ES.C_cur_load_gas,
        ES.π, ES.compressors, ES.C, ES.S, ES.D_gas, ES.T)
    ES.gas_cost = total_cost_gas(
        ES.objective_type, ES.dt, ES.q_s, ES.sources, ES.q_d_cur, ES.C_cur_load_gas,
        ES.π, ES.compressors, ES.C, ES.S, ES.D_gas, ES.T)
end


function get_flow_directions(m::Matrix{Float64})
    # if model_type == "trade"
    z = ones(size(m))
    z[m .< 0] .= 0
    # end
    return z
end


function calc_γ(m, π_avg)
    return m*abs(m)/π_avg
end


function write_solution_power!(
    IES::IntegratedEnergySystem, model::JuMP.Model, model_type::String)
    """
    Extracts and writes model solutions for the power system side on the
    IntegratedEnergySystem object.
    """

    write_solution_vars!(IES, model, [:p, :w, :θ, :p_d_cur])
    GFPPs = [g for g in IES.G if isa(IES.generators[g], GasFiredGenerator)]
    IES.q_g = [IES.generators[g].Γ*value(model[:p][g,t]) for g in GFPPs, t in IES.T]
    IES.f = [1/IES.lines[l].X*(IES.θ[IES.lines[l].start,t]-IES.θ[IES.lines[l].stop,t])
        for l in IES.L, t in IES.T]

    nonGFPPs = [g for g in IES.G if isa(IES.generators[g], NonGasFiredGenerator)]
    IES.power_cost = cost_power_linear(
        IES.dt, IES.p, IES.generators, nonGFPPs, IES.p_d_cur, IES.D_el, IES.C_cur_load_power, IES.T)
    IES.total_cost += IES.power_cost
end


function write_solution!(GS::GasSystem, model::JuMP.Model, model_type::String)
    write_solution_gas!(GS, model, model_type)
    return nothing
end


function write_solution!(IES::IntegratedEnergySystem, model::JuMP.Model, model_type::String)
    write_solution_gas!(IES, model, model_type)
    write_solution_power!(IES, model, model_type)
    return nothing
end
