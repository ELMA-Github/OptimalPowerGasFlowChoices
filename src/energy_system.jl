abstract type EnergySystem end

mutable struct GasSystem <: EnergySystem
    path_input_data_gas :: String # input folder of scenario data
    case_study :: String # Case study
    nodes #:: Dict{Int64, Node}() # Node container
    sources #:: Dict{Int64, Source}() # Source container
    demands #:: Dict{Int64, Any}() # Demand container -- can be either type Gasload_uniform or Gasload_average
    pipes #:: Dict{Int64, Pipeline}() # Pipe container
    compressors #:: Dict{Int64, Compressor}() # Compressor container
    I # Set of nodes
    S # Set of sources
    P # Set of pipes
    P_hat # Set of undiscretized pipes in the network
    D_gas # Set of gas demands
    C # Set of compressors
    space_disc :: Int64 # Control parameter if space discretization should be performed
    segment :: Int64 # Length of each subpipeline segment
    timehorizon :: Int64 # Time horizon (hours)
    print_model :: Int64 # Controls if model is printed
    T # Set of timesteps
    dt :: Int64 # Seconds per timestep
    N_dt :: Int64 # Number of timesteps in the time horizon#
    PDE_type :: Int # Control parameter for PDE model
    model_type :: String # Model type
    objective_type :: String # Linear or quadratic cost function
    C_cur_load_gas :: Number # Value of lost load (\$/kgh)
    obj_val # Objective function value of model
    total_cost # Optimal total cost
    gas_cost # Optimal total gas cost
    solve_time # Solve time of model instance
    q_s # Optimal gas production
    π # Optimal nodal pressure
    m # Optimal gas flow
    m_in # Optimal gas flow in
    m_out # Optimal gas flow Q_nm_out
    π_avg # Optimal average pressure in pipelines
    q_d_cur # Optimal curtailment
    m_c # Optimal gas flow through compressors
    LP # Optimal linepack mass
    z # Optimal flow directions
    m⁺ # Optimal average flow in pipeline nm in postive direction
    m⁻ # Optimal average flow in pipeline nm in postive direction
    γ
    γ⁺
    γ⁻
    pr_restored # restored pressure values from original PDE problem
    obj_val_restored # restored objective function value from original PDE problem
    config_dict_solver # Tuple, contains (1) solver and (2) solver attribute dictionary

    model # Optimization model
    warmstart::Dict{Symbol, Any}

    pu::Bool
    bases_dict::Dict{Symbol,Number}
    boundary_conditions::Dict{Symbol,Vector}

    # Adapted from https://discourse.julialang.org/t/default-value-of-some-fields-in-a-mutable-struct/33408/20
    # Authors: GunnarFarneback, davidbp
    function GasSystem(timehorizon, dt, case_study; kwargs...)
        GS = new()
        GS.N_dt = Int64(timehorizon*3600/dt)
        GS.T = sort(collect(1:GS.N_dt))
        GS.path_input_data_gas = joinpath(dirname(@__DIR__), "inputs", "gas_only", case_study)
        for (key, value) in kwargs
            # field_type_key = typeof(getfield(GS, key))
            setfield!(GS, key, value)
            # setfield!(GS, key, convert(field_type_key, value))
        end
        return GS
    end
end


function get_gas_system_attributes(ES::EnergySystem)

    return ES.nodes, ES.sources, ES.demands, ES.pipes, ES.compressors,
        ES.I, ES.S, ES.D_gas, ES.P, ES.P_hat, ES.T, ES.C, ES.boundary_conditions,
        ES.dt, ES.C_cur_load_gas
end


mutable struct IntegratedEnergySystem <: EnergySystem
    ##### Gas system attributes #####
    path_input_data_gas :: String # input folder of gas scenario data
    case_study :: String # Case study
    nodes #:: Dict{Int64, Node}() # Node container
    sources #:: Dict{Int64, Source}() # Source container
    demands #:: Dict{Int64, Any}() # Demand container -- can be either type Gasload_uniform or Gasload_average
    pipes #:: Dict{Int64, Pipeline}() # Pipe container
    compressors #:: Dict{Int64, Compressor}() # Compressor container
    I # Set of nodes
    S # Set of sources
    P # Set of pipes
    P_hat # Set of undiscretized pipes in the network
    D_gas # Set of gas demands
    C # Set of compressors
    space_disc :: Int64 # Control parameter if space discretization should be performed
    segment :: Int64 # Length of each subpipeline segment
    timehorizon :: Int64 # Time horizon (hours)
    print_model :: Int64 # Controls if model is printed
    T # Set of timesteps
    dt :: Int64 # Seconds per timestep
    N_dt :: Int64 # Number of timesteps in the time horizon#
    PDE_type :: Int # Control parameter for PDE model
    model_type :: String # Model type
    objective_type :: String # Linear or quadratic cost function
    C_cur_load_gas :: Number # Value of lost load (\$/kgh)
    obj_val # Objective function value of model
    total_cost # Optimal total cost
    gas_cost # Optimal total gas cost
    power_cost # Optimal total power cost
    solve_time # Solve time of model instance
    q_s # Optimal gas production
    π # Optimal nodal pressure
    m # Optimal gas flow
    m_in # Optimal gas flow in
    m_out # Optimal gas flow Q_nm_out
    π_avg # Optimal average pressure in pipelines
    q_d_cur # Optimal curtailment
    m_c # Optimal gas flow through compressors
    LP # Optimal linepack max_optimal_physics_gap
    z # Optimal flow directions
    m⁺ # Optimal average flow in pipeline nm in postive direction
    m⁻ # Optimal average flow in pipeline nm in postive direction
    γ
    γ⁺
    γ⁻
    pr_restored # restored pressure values from original PDE problem
    obj_val_restored # restored objective function value from original PDE problem

    ##### Power system attributes #####
    path_input_data_power :: String # input folder of gas scenario data
    buses # Dict containing power bus objects
    lines # Dict containing power line objects
    generators # Dict containing dispatchable generator objects
    windgenerators # Dict containing wind power objects
    demands_EL # Dict containing power demand objects
    N # List of buses
    L # List of lines
    G # List of dispatchable generators
    J # List of wind power producers
    D_el # List of power demands
    q_g # Gas consumption of GFPPs at optimum
    p # Generator dispatch at optimum
    w # Wind generator dispatch at optimum
    θ # Bus voltage angles at optimum
    f # Line power flow at optimum
    p_d_cur # Power curtailment at optimum
    C_cur_load_power :: Number # Value of lost load (\$/MW)
    config_dict_solver # Tuple, contains (1) solver and (2) solver attribute dictionary

    model # Optimization model
    warmstart::Dict{Symbol, Any}

    pu::Bool
    bases_dict::Dict{Symbol,Number}
    boundary_conditions::Dict{Symbol,Vector}

    # Adapted from https://discourse.julialang.org/t/default-value-of-some-fields-in-a-mutable-struct/33408/20
    # Authors: GunnarFarneback, davidbp
    function IntegratedEnergySystem(timehorizon, dt, case_study; kwargs...)
        IES = new()
        IES.N_dt = Int64(timehorizon*3600/dt)
        IES.T = sort(collect(1:IES.N_dt))
        IES.path_input_data_gas = joinpath(dirname(@__DIR__), "inputs", "integrated_power_and_gas", case_study, "gas")
        IES.path_input_data_power = joinpath(dirname(@__DIR__), "inputs", "integrated_power_and_gas", case_study, "power")
        for (key, value) in kwargs
            # field_type_key = typeof(getfield(ES, key))
            setfield!(IES, key, value)
            # setfield!(ES, key, convert(field_type_key, value))
        end
        return IES
    end
end

function get_power_system_attributes(ES::IntegratedEnergySystem)

    return ES.generators, ES.windgenerators, ES.demands_EL, ES.lines,
        ES.N, ES.G, ES.J, ES.G, ES.D_el, ES.L, ES.buses, ES.C_cur_load_power
end


function intialize_energy_system(
    config_dict::Dict{Symbol,Any}
    )

    system = pop!(config_dict, :system)
    if system == "gas_only"
        ES = GasSystem(
            config_dict[:timehorizon], 
            config_dict[:dt],
            config_dict[:case_study];
            config_dict...)

    elseif system == "integrated_power_and_gas"
        ES = IntegratedEnergySystem(
            config_dict[:timehorizon], 
            config_dict[:dt],
            config_dict[:case_study];
            config_dict...)
    else 
        throw(ArgumentError("System is set to $(config_dict[:system]).
            Must be either 'integrated_power_and_gas' or 'gas_only'."))
    end
    return ES
end

fieldnames_gassystem = [
    :case_study, :nodes, :sources, :demands,
    :pipes, :compressors, :I, :S, :P, :D_gas, :C, :space_disc, :segment,
    :timehorizon, :print_model, :dt, :PDE_type,
    :model_type, :C_cur_load_gas
    ] 
fieldnames_integratedenergysystem = vcat(fieldnames_gassystem,
    [:buses, :lines, :generators, :windgenerators,
    :demands_EL, :N, :L, :G, :J, :D_el])
    
Base.copy(GS::GasSystem) = GasSystem(GS.timehorizon, GS.dt, GS.case_study;
    Dict(k => getfield(GS, k) for k in fieldnames_gassystem)...)

Base.copy(IES::IntegratedEnergySystem) = IntegratedEnergySystem(
    IES.timehorizon, IES.dt, IES.case_study;
    Dict(k => getfield(IES, k) for k in fieldnames_integratedenergysystem)...)