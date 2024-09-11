include("energy_system.jl")
include("PTDF_matrix.jl")
include("bases_choice.jl")


##================= Gas base values ===============
function create_base_gas(path_input_data::String, pu::Bool=true, dt::Number=3600, segment::Number=50000, model_type::AbstractString="NCNLP")

    @info "Create pu base."
    if pu
        Base_df = CSV.read(joinpath(
            path_input_data, "gas_params.csv"), DataFrame, delim=",")
        Π_base, D_base, T_base, Cost_base = Base_df[1, [
            :Pr_base_MPa, :D_base_m, :T_base_s, :Cost_base_dollar]]
        @info "User choice:"
    else
        Π_base, D_base, T_base, Cost_base = 1, 1e-6, 1, 1000
        @info "Default choice (no pu):"
    end
    # @info "Π_base: $Π_base, D_base: $D_base, T_base: $T_base, Cost_base: $Cost_base"
    @info "Π_base: $Π_base, Cost_base: $Cost_base"

    # c = model_type in ["boundary_conditions", "NCNLP"] ? -1 : 0
    c = 1
    # D_base, T_base = find_optimal_bases(0.5, dt, segment, c)
    # pipe_df = CSV.read(joinpath(path_input_data, "gas_pipes.csv"), DataFrame, delim=",")
    @info "D_base: $D_base, T_base: $T_base"
    @info "Calculating residual bases."

    # Residual bases
    A_base = D_base^2 # [m^2]
    # ρ_base = Π_base * 10^6 / C_base^2 # [kg/m^3]
    # U_base = M_base/(ρ_base*A_base) # [m/s]

    # L_base and T_base can be derived by conservation of mass and momentum constraints
    Π_base_Pa = Π_base*10^6
    C_base = D_base/T_base
    M_base = D_base*Π_base_Pa*T_base
    Γ_base = Π_base_Pa*D_base^2*T_base^2
    K_base = D_base*T_base

    @info "Residual bases:"
    @info "A_base: $A_base, C_base: $C_base, M_base: $M_base, Γ_base: $Γ_base, K_base: $K_base"

    bases_dict = Dict{Symbol, Number}(
        :Π_base => Π_base, # [MPa]=[Mkg/(m*s^2)]
        :M_base => M_base, # [kg/s]
        :D_base => D_base, # [m]
        :C_base => C_base, # [m/s]
        :T_base => T_base, # [s]
        :K_base => K_base, # [m^2*s^2]
        :Cost_base => Cost_base, # [$]
        :Γ_base => Γ_base # [kg*m]
    )

    return bases_dict
end


##================= Power base values ===============
function create_base_power(path_input_data::String, pu::Bool=true)

    if pu
        Base_df = CSV.read(joinpath(
            path_input_data, "el_params.csv"), DataFrame, delim=",")
        S_base = Base_df[1, [
            :S_base_MVA]][1]
    else
        S_base = 100
    end

    bases_dict = Dict{Symbol, Number}(
        :S_base => S_base, # [MVA]
    )

    return bases_dict
end


function unpack_gas_bases_dict(bases_dict::Dict{Symbol,Number})
    bases_list = [:Π_base, :M_base, :D_base, :C_base, :Cost_base]
    return [bases_dict[bases] for bases in bases_list]
end


function unpack_power_bases_dict(bases_dict::Dict{Symbol,Number})
    bases_list = [:S_base]
    return [bases_dict[bases] for bases in bases_list]
end


function unpack_bases_dict(bases_dict::Dict{Symbol,Number})
    return vcat(unpack_gas_bases_dict(bases_dict), unpack_power_bases_dict(bases_dict))
end


##================= Gas nodes data ===============
function create_Nodes(path_input_data::String, Π_base::Number)
    N_df = CSV.read(joinpath(path_input_data, "gas_nodes.csv"), DataFrame, delim=",")
    nodes = Dict{Int,Node}()
    for row in 1:nrow(N_df)

        id, Π_max, Π_min, node_type, Π_slack =
            N_df[row, [:Node_No, :Pmax_MPa, :Pmin_MPa, :Node_Type, :Pslack_MPa]]

        x, y = 0, 0
        try 
            x, y = N_df[row, [:x, :y]]
        catch error
        end

        if node_type == 0
            nodes[id] =  NonSlackNode(id, Π_max, Π_min, x, y, Π_base) 
        elseif node_type == 1
            nodes[id] =  SlackNode(id, Π_max, Π_min, Π_slack, x, y, Π_base) 
        else throw(ArgumentError("Node_type is $node_type. Must be either '1' or '0'."))
        end
    end
    return nodes
end


##================= Gas sources data ===============
function create_GasSources(path_input_data::String, M_base::Number, Cost_base::Number)
    S_df = CSV.read(joinpath(path_input_data, "gas_supply.csv"), DataFrame, delim=",")
    sources = Dict{Int,GasSource}()
    for row in 1:nrow(S_df)
        id, location, Q_max, Q_min, C_1, C_2 =
            S_df[row, [:Supply_No, :Node, :Smax_kg_s, :Smin_kg_s,:C1_per_kgh, :C2_per_kgh2]]
        sources[id] =  GasSource(id, location, Q_max, Q_min, C_1, C_2, M_base, Cost_base) 
    end
    return sources
end



function isuniform(
    dt_data::Int, # Resolution of input data [s]
    dt_disc::Int # Time interval discretization [s])
    )
    """
    Determines whether the load is uniformly distributed (true) or averaged (false)
    over a time interval of the model depending on the given data.
    """
    return dt_data >= dt_disc 
end

##================= Gas load data===============
function create_GasDemands(
    path_input_data::String, dt::Int, dt_gasload_s::Int, M_base::Number)
    D_df = CSV.read(joinpath(path_input_data, "gas_load.csv"), DataFrame, delim=",")
    demands = Dict{Int,GasLoad}()
    GP_df= CSV.read(joinpath(path_input_data, "gas_profile.csv"), DataFrame, delim=",")
    unif = isuniform(dt_gasload_s, dt)
    for row in 1:nrow(D_df)
        id, node_load, load_nom, type_profile = 
            D_df[row, [:Load_No,:Node,:Load_kg_s,:Profile]]
        load_profile = GP_df[:,type_profile]
        demands[id] =  GasLoad(
            id, node_load, load_nom, type_profile, load_profile, dt, dt_gasload_s, unif, M_base) 
    end
    return demands
end


##================= Gas pipeline data ===============
function create_Pipes(
    path_input_data::String, nodes::Dict{Int,Node}, D_base::Number,
    C_base::Number)
    P_df = CSV.read(joinpath(path_input_data, "gas_pipes.csv"), DataFrame, delim=",")
    pipes = Dict{Int,Pipeline}()
    for row in 1:nrow(P_df)
        id, start, stop, F, D, L =
            P_df[row, [:Pipe_No, :From_Node, :To_Node, :friction, :Diameter_m, :Length_m]]
         pipes[id] =  Pipeline(
            id, start, stop, id, F, D, L, nodes, D_base, C_base)
    end
    return pipes
end


##================= Gas compressor data ===============
function create_Compressors(
    path_input_data::String,
    # nodes::Dict{Int,Node},
    Π_base::Number,
    Cost_base::Number)
    C_df = CSV.read(joinpath(path_input_data, "gas_compressors.csv"), DataFrame, delim=",")
    compressors = Dict{Int,Compressor}()
    for row in 1:nrow(C_df)
        id, from_node, to_node, fuel_node, fuel_consumption, CR_max, CR_min, compr_cost =
            C_df[row, [:Compressor_No,:From_Node,:To_Node,:fuel_gas_node,
                :fuel_gas_consumption,:CR_Max,:CR_Min,:Compression_cost]]
        compressors[id] =  Compressor(
            id, from_node, to_node, fuel_node, fuel_consumption, CR_max, CR_min, compr_cost, 
            # nodes,
            Π_base, Cost_base) 
    end
    return compressors
end


##================= Electricity bus data ===============
function create_Buses(path_input_data)
    B_df = CSV.read(joinpath(path_input_data, "buses_EL.csv"), DataFrame, delim=",")

    buses = Dict{Int,Bus}()
    for row in 1:nrow(B_df)
        id, slack  = B_df[row, [:Bus_No, :Slack]]
        buses[id] =  Bus(id, slack) 
    end
    return buses
end


##================= Electricity DispatchableGenerators data ===============
function create_DispatchableGenerators(
    path_input_data::AbstractString, S_base::Number, M_base::Number, Cost_base::Number)
    G_df = CSV.read(joinpath(path_input_data, "dispatchablegenerators.csv"), DataFrame, delim=",")
    generators = Dict{Int,DispatchableGenerator}()
    for row in 1:nrow(G_df)
        id, bus, P_min, P_max, P_down, P_up, type, node, eta, c1_el, c2_el = 
            G_df[row, [:Gen_num,:EL_node,:Pmin_MW,:Pmax_MW,:P_down_MW_h,
            :P_up_MW_h,:Type,:NG_node,:Conversion_kg_sMW,:C1_per_MWh,:C2_per_MWh2]]
        if type == "non-NGFPP"
            generators[id] =  NonGasFiredGenerator(
                id, bus, P_min, P_max, P_down, P_up, c1_el, c2_el, S_base, Cost_base) 
        elseif type == "NGFPP"
            generators[id] =  GasFiredGenerator(
                id, bus, P_min, P_max, P_down, P_up, Int64(node), eta, S_base, M_base)
        else
            throw(AttributeError("""Type if dispatchable power plant is $type.
            Only types 'non-NGFPP' and 'NGFPP' are allowed."""))
        end
    end
    return generators
end


function check_validty_timehorizon(horizon_model::Int, horizon_data::Int)
    if horizon_model > horizon_data
        throw(ArgumentError("""
            Input timehorizon $horizon_model hours is longer than the available data.
            Timehorizons larger than $horizon_data are not allowed for this case study."""))
    end
    return nothing
end


function ismultiple(val1, val2)
    return val1 % val2 == 0
end


function check_validity_time_discretization(dt_disc::Int, dt_data::Int, var::AbstractString)
    if dt_disc >= dt_data
        if ~ismultiple(dt_disc, dt_data)
            throw(ArgumentError("""
                Time discretization is $dt_disc while resolution of $var data is $dt_data.
                If the time discretization shall be larger than the resolution of the data, it must be a whole multiple.
                For example, $(2*dt_data) or $(4*dt_data)."""))
        end
    else
        if ~ismultiple(dt_data, dt_disc)
            throw(ArgumentError("""
                Time discretization is $dt_disc while resolution of $var data is $dt_data.
                If the time discretization shall be smaller than the resolution of the data, it must be a whole multiple.
                For example, $(Int64(0.5*dt_data)) or $(Int64(0.25*dt_data))."""))
        end
    end
    return nothing
end


##================= Gas base values ===============
function create_gas_input_params(path_input_data::String, timehorizon::Int, dt_disc::Int)

    Base_df = CSV.read(joinpath(
        path_input_data, "gas_params.csv"), DataFrame, delim=",")

    T_gasload_h, dt_gasload_s = Base_df[1, [
        :T_gasload_h, :dt_gasload_s]]

    check_validty_timehorizon(timehorizon, T_gasload_h)
    check_validity_time_discretization(dt_disc, dt_gasload_s, "gas load")

    return dt_gasload_s
end


##================= Electricity base values ===============
function create_el_input_params(path_input_data::String, timehorizon::Int, dt_disc::Int)

    Base_df = CSV.read(joinpath(
        path_input_data, "el_params.csv"), DataFrame, delim=",")

    T_eload_h, dt_eload_s, T_wind_h, dt_wind_s = Base_df[1, [
        :T_eload_h, :dt_eload_s, :T_wind_h, :dt_wind_s]]

    T_max = min(T_eload_h, T_wind_h) 
    check_validty_timehorizon(timehorizon, T_max)
    for dt_data in [(dt_eload_s, "electricity load"), (dt_wind_s, "wind generation")]
        check_validity_time_discretization(dt_disc, dt_data[1], dt_data[2])
    end

    return dt_eload_s, dt_wind_s
end

##================= Electricity WindGenerators data ===============
function create_WindGenerators(path_input_data::AbstractString, dt::Number, dt_wind_s, S_base::Number)
    W_df = CSV.read(joinpath(path_input_data, "windgenerators.csv"), DataFrame, delim=",")
    windgenerators = Dict{Int,WindGenerator}()
    WP_df = CSV.read(joinpath(path_input_data, "wind_profile.csv"), DataFrame, delim=",")
    unif = isuniform(dt_wind_s, dt)
    for row in 1:nrow(W_df)
        id, bus, P_max, type_profile = W_df[row, [:Wind_num,:EL_node,:Pmax_MW,:profile_type]]
        wind_profile = WP_df[:,type_profile]
        windgenerators[id] =  WindGenerator(
            id, bus, P_max, type_profile, wind_profile, dt, dt_wind_s, unif, S_base) 
    end
    return windgenerators
end

##================= Electricity demand data ===============
function create_PowerDemands(path_input_data::AbstractString, dt::Number, dt_eload_s, S_base::Number)
    D_el_df = CSV.read(joinpath(path_input_data, "electricity_load.csv"), DataFrame, delim=",")
    demand_EL = Dict{Int,PowerDemand}()
    EP_df = CSV.read(joinpath(path_input_data, "electricity_profile.csv"), DataFrame, delim=",")
    unif = isuniform(dt_eload_s, dt)
    for row in 1:nrow(D_el_df)
        id, bus, load_nom, type_profile = D_el_df[row,
            [:Load_No,:EL_Node,:Load_MW,:Profile]]
        load_profile = EP_df[:,type_profile]
        demand_EL[id] =  PowerDemand(id, bus, load_nom, load_profile, dt, dt_eload_s, unif, S_base)  
    end
    return demand_EL
end


##================= Electricity Line data ===============
function create_Lines(path_input_data::AbstractString, S_base::Number)
    L_df = CSV.read(joinpath(path_input_data, "lines.csv"), DataFrame, delim=",")
    Ψ = create_PTDF_matrix(L_df)
    Ψ = round.(Ψ, digits=4)

    lines= Dict{Int,Line}()
    for row in 1:nrow(L_df)
        id, start, stop, X, F_max = L_df[row, [:Line_num,:Start,:Stop,:X_pu,:Capacity_MW]]
        lines[id] =  Line(id, start, stop, X, F_max, Ψ[row,:], S_base) 
    end
    return lines
end


function get_ordered_sets(components::Vector{Dict{Int}})
    """
    Creates ordered sets of each of the passed components,
    e.g., Pipelines, which are stored in a pipes dictionray.
    """
    return [sort(collect(keys(c))) for c in components]
end


function initialize_gas_system_data(
    path_input_data_gas, dt, timehorizon, bases_dict, space_disc, segment, pu)

    Π_base, M_base, D_base, C_base, Cost_base =
        unpack_gas_bases_dict(bases_dict)

    dt_gasload_s = create_gas_input_params(
        path_input_data_gas, timehorizon, dt)

    ##================= Instantiate gas network data ===============
    nodes = create_Nodes(path_input_data_gas, Π_base)
    pipes = create_Pipes(path_input_data_gas, nodes, D_base, C_base)
    sources = create_GasSources(path_input_data_gas, M_base, Cost_base) 
    demands = create_GasDemands(path_input_data_gas, dt, dt_gasload_s, M_base)
    compressors = create_Compressors(path_input_data_gas, Π_base, Cost_base)

    ##================= Run space discretization gas network ===============
    if space_disc == 1
        sub_nodes, sub_pipes = space_discretization(nodes, pipes, segment, bases_dict)
    else
        sub_nodes, sub_pipes = nodes, pipes
    end

    P_hat = sort(collect(keys(pipes)))

    # read_boundary_conditions(path_input_data_gas)
    boundary_conditions = Dict()

    return sub_nodes, sub_pipes, sources, demands, compressors, boundary_conditions, P_hat
end


function initialize_power_system_data(
    path_input_data_power, dt, timehorizon, bases_dict)

    S_base, M_base, Cost_base = 
        bases_dict[:S_base], bases_dict[:M_base], bases_dict[:Cost_base]

    dt_eload_s, dt_wind_s = create_el_input_params(
        path_input_data_power, timehorizon, dt)

    ##================= Instantiate power network data ===============
    buses = create_Buses(path_input_data_power)
    lines = create_Lines(path_input_data_power, S_base)
    generators = create_DispatchableGenerators(
        path_input_data_power, S_base, M_base, Cost_base)
    windgenerators = create_WindGenerators(path_input_data_power, dt, dt_wind_s, S_base)
    demands_EL = create_PowerDemands(path_input_data_power, dt, dt_eload_s, S_base)
    
    return buses, lines, generators, windgenerators, demands_EL
end


function intialize_data!(GS::GasSystem)
    """
    Loads and returns input data from CSV files for a given case study.
    """

    ##================= Instantiate bases ===============
    GS.bases_dict = create_base_gas(GS.path_input_data_gas, GS.pu, GS.dt, GS.segment, GS.model_type)
    GS.C_cur_load_gas = GS.C_cur_load_gas*GS.bases_dict[:M_base]/GS.bases_dict[:Cost_base] 

    ##================= Instantiate gas system data ===============
    GS.nodes, GS.pipes, GS.sources, GS.demands, GS.compressors, GS.boundary_conditions, GS.P_hat =
        initialize_gas_system_data(
            GS.path_input_data_gas, GS.dt, GS.timehorizon,
            GS.bases_dict, GS.space_disc, GS.segment, GS.pu)

    ##================= Create ordered sets ===============
    GS.I, GS.S, GS.P, GS.D_gas, GS.C = get_ordered_sets(
        [GS.nodes, GS.sources, GS.pipes, GS.demands, GS.compressors])

    return nothing

end


function intialize_data!(IES::IntegratedEnergySystem)
    """
    Loads and returns input data from CSV files for a given case study.
    """

    ##================= Instantiate bases ===============
    bases_dict_gas = create_base_gas(IES.path_input_data_gas, IES.pu, IES.dt, IES.segment, IES.model_type)
    bases_dict_power = create_base_power(IES.path_input_data_power, IES.pu)
    IES.bases_dict = merge(bases_dict_gas, bases_dict_power)
    IES.C_cur_load_gas = IES.C_cur_load_gas*IES.bases_dict[:M_base]/IES.bases_dict[:Cost_base] 
    IES.C_cur_load_power = IES.C_cur_load_power*IES.bases_dict[:S_base]/IES.bases_dict[:Cost_base] 

    ##================= Instantiate gas system data ===============
    IES.nodes, IES.pipes, IES.sources, IES.demands, IES.compressors, IES.boundary_conditions, IES.P_hat =
        initialize_gas_system_data(
            IES.path_input_data_gas, IES.dt, IES.timehorizon,
            IES.bases_dict, IES.space_disc, IES.segment, IES.pu)

    # ##================= Instantiate power system data ===============
    IES.buses, IES.lines, IES.generators, IES.windgenerators, IES.demands_EL =
        initialize_power_system_data(IES.path_input_data_power, IES.dt, IES.timehorizon, IES.bases_dict)

    ##================= Create ordered sets ===============
    ##### Gas #####
    IES.I, IES.S, IES.P, IES.D_gas, IES.C = get_ordered_sets(
        [IES.nodes, IES.sources, IES.pipes, IES.demands, IES.compressors])

    ##### Power #####
    IES.N, IES.L, IES.G, IES.J, IES.D_el = get_ordered_sets(
        [IES.buses, IES.lines, IES.generators, IES.windgenerators, IES.demands_EL])

    return nothing
end


function PDE_type_parameters(PDE_type::Int)
    """
    Chooses model parameters ut1 and ut2 for the PDEs of the gas flow model.

    Returns:
    Tuple, First entry is ut1, second entry ut2
    """

    PDE_type_dict = Dict(
        1 => (0,0), #steady-state model
        2 => (1,0), #quasi-dynamic model
        3 => (1,1)  #transient model
    )

    return PDE_type_dict[PDE_type]
end