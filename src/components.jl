############### Gas network structs ###############
struct GasSource
    id::Int # id
    node::Int # Gas source location
    Q_s_max::Number # Maximum source [p.u.]
    Q_s_min::Number # Minimum source [p.u.]
    C_1::Number # Linear cost coefficient [p.u.]
    C_2::Number # Quadratic cost coefficient [p.u.]
    function GasSource(
        id::Int, # id
        node::Int, # Gas source location
        Q_s_max::Number, # Maximum source [kg/s]
        Q_s_min::Number, # Minimum source [kg/s]
        C_1::Number, # Linear cost coefficient [$/kg]
        C_2::Number, # Quadratic cost coefficient [$/(kg)^2]
        M_base::Number, # Mass flow base [kg/s]
        Cost_base::Number, # Cost base [$]
    )
        Q_s_max, Q_s_min = Q_s_max/M_base, Q_s_min/M_base
        C_1, C_2 = C_1*M_base/Cost_base, C_2*M_base^2/Cost_base

        return new(id, node, Q_s_max, Q_s_min, C_1, C_2)
    end
end


abstract type Node end

struct NonSlackNode <: Node
    id::Int # id
    Π_max::Number # Maximum nodal pressure [p.u.]
    Π_min::Number # Minimum nodal pressure [p.u.]
    x::Number # x-Coordinate for illustrations
    y::Number # y-Coordinate for illustrations
    function NonSlackNode(
        id::Int, # id
        Π_max::Number, # Maximum nodal pressure [MPa]
        Π_min::Number, # Minimum nodal pressure [MPa]
        x::Number, # x-Coordinate for illustrations
        y::Number, # y-Coordinate for illustrations
        Π_base::Number=1, # Pressure base [MPa]
    )
        
        Π_max, Π_min = Π_max/Π_base, Π_min/Π_base
        return new(id, Π_max, Π_min, x, y)
    end
end


struct SlackNode <: Node
    id::Int # id
    Π_max::Number # Maximum nodal pressure [p.u.]
    Π_min::Number # Minimum nodal pressure [p.u.]
    Π_slack::Number # 1:slack node, 0:otherwise
    x::Number # x-Coordinate for illustrations
    y::Number # y-Coordinate for illustrations
    function SlackNode(
        id::Int, # id
        Π_max::Number, # Maximum nodal pressure [MPa]
        Π_min::Number, # Minimum nodal pressure [MPa]
        Π_slack::Number, # Slack pressure [MPa]
        x::Number, # x-Coordinate for illustrations
        y::Number, # y-Coordinate for illustrations
        Π_base::Number=1, # Pressure base [MPa]
    )
        
        Π_max, Π_min, Π_slack = Π_max/Π_base, Π_min/Π_base, Π_slack/Π_base
        return new(id, Π_max, Π_min, Π_slack, x, y)
    end
end


const SPEED_OF_SOUND = 350 # Speed of sound in gas [m/s]


struct Pipeline
    id::Int # id
    start::Int # Start node id
    stop::Int # End node id
    parent::Int # id of parent pipe
    F::Number # Friction coefficient [dimensionless]
    C::Number # Speed of sound [p.u.]
    D::Number # Diameter [p.u.]
    L::Number # Length [p.u.]
    A::Number # Cross sectional area [p.u.]
    S::Number # Linepack parameter [?]
    LP0::Number # initial line pack mass [?]
    K::Number # Weymouth coefficient [p.u.]
    M_min::Number # Min average flow rate in pipelines [p.u.]
    M_max::Number # Max average flow rate in pipelines [p.u.]
    Π_avg_min::Number # Minimum average nodal pressure [p.u.]
    Π_avg_max::Number # Maximum average nodal pressure [p.u.]
    Γ_min::Number # Minimum bound on momentum conservation auxiliary [p.u.]
    Γ_max::Number # Maximum bound on momentum conservation auxiliary [p.u.]
    x::Number # x-Coordinate for illustrations (pipeline midpoint)
    y::Number # y-Coordinate for illustrations (pipeline midpoint)
    function Pipeline(
        id::Int, # id
        start::Int, # Start node id
        stop::Int, # End node id
        parent::Int, # id of parent pipe
        F::Number, # Friction coefficient [dimensionless]
        D::Number, # Diameter [m]
        L::Number, # Length [m]
        nodes::Dict{Int,<:Node}, # Nodes dict
        D_base::Number, # Diameter base [m]
        C_base::Number, # Speed of sound in gas base [m/s]
    )

        D = D / D_base
        C = SPEED_OF_SOUND / C_base
        L = L / D_base

        A = pi*(D/2)^2
        K = sqrt(D*A^2/(F*C^2*L))
        S = L*A/C^2
        LP0 = L*A/C^2 *1 # 1 is arbitrary for average pressure

        # Π_min and Π_max are in pu
        Π_avg_min = (nodes[start].Π_min+nodes[stop].Π_min)/2
        Π_avg_max = (nodes[start].Π_max+nodes[stop].Π_max)/2
        M_min = -K*sqrt(nodes[stop].Π_max^2-nodes[start].Π_min^2) # Approximate by Weymouth-Equation
        M_max = K*sqrt(nodes[start].Π_max^2-nodes[stop].Π_min^2) # Approximate by Weymouth-Equation

        Γ_min = M_min*abs(M_min)/((nodes[stop].Π_max+nodes[start].Π_min)/2)
        Γ_max = M_max*abs(M_max)/((nodes[start].Π_max+nodes[stop].Π_min)/2)

        x = (nodes[start].x+nodes[stop].x)/2
        y = (nodes[start].y+nodes[stop].y)/2

        return new(id, start, stop, parent, F, C, D, L, A, S, LP0, K,
            M_min, M_max, Π_avg_min, Π_avg_max, Γ_min, Γ_max, x, y)
    end
end


struct Compressor
    id::Int # id
    start::Int # Start node id
    stop::Int # End node id
    fuel_node::Int # Fuel gas node id
    fuel_consumption::Number # Fuel gas consumption [dimensionless]
    CR_max::Number # Maximum compression ratio [dimensionless]
    CR_min::Number # Minimum compression ratio [dimensionless]
    C_compr::Number # Cost compression [p.u.]
    function Compressor(
        id::Int, # id
        start::Int, # Start node id
        stop::Int, # End node id
        fuel_node::Int, # Fuel gas node id
        fuel_consumption::Number, # Fuel gas consumption [dimensionless]
        CR_max::Number, # Maximum compression ratio [dimensionless]
        CR_min::Number, # Minimum compression ratio [dimensionless]
        C_compr::Number, # # Cost compression [cost/MPa]
        Π_base::Number, # Pressure base [MPa]
        Cost_base::Number # Cost base [$]
    )
    """
    Gas compressor that consumes gas as a fuel to power its drive from fuel_node.
    The fuel consumption is determined by the fraction fuel_consumption of the compressed gas flow.
    """

        C_compr = C_compr/Cost_base*Π_base
        return new(
            id, start, stop, fuel_node, fuel_consumption, CR_max, CR_min, C_compr)
    end
end


struct GasLoad
    id::Int # id
    node::Int # Gas node id
    type_load::String # Gas load type {uniform, average}
    load_nom::Number # Nominal load [p.u.]
    type_profile::AbstractString # Type of gas load profile
    load_profile::Vector{Number} # Gas load profile (-) [dimensionless]
    Q_d::Vector{Number} # Gas load at each timestep [p.u.]
    function GasLoad(
        id::Int, # id
        node::Int, # Gas node id
        load_nom::Number, # Nominal load [kg/s]
        type_profile::AbstractString, # Type of gas load profile
        load_profile::Vector{<:Number}, # Gas load profile (-) [dimensionless]
        dt::Int, # Time interval discretization [s]
        dt_gasload_s::Int, # Resolution of input data [s]
        uniform::Bool, # Uniform (true) or average (false) profile
        M_base::Number, # Mass flow base [kg/s]
    )

        load_nom = load_nom/M_base
        type, Q_d = create_profile(uniform, load_nom, load_profile, dt_gasload_s, dt)

        return new(id, node, type, load_nom, type_profile, load_profile, Q_d)
    end
end


function uniform_profile(
    nom_val::Number, # Nominal value [p.u.]
    profile::Vector{<:Number}, # Nominalized profile
    dt_data::Int, # Resolution of input data [s]
    dt_disc::Int # Time interval discretization [s]
    )
    """
    Assumes a uniformly distributed profile when expanding a time series.
    Thus, the profile does not change between steps of the given time interval of the model.
    """
    val = repeat(nom_val.*profile, inner=Int64(dt_data/dt_disc))
    return val
end


function average_profile(
    nom_val::Number, # Nominal value [p.u.]
    profile::Vector{<:Number}, # Nominalized profile
    dt_data::Int, # Resolution of input data [s]
    dt_disc::Int # Time interval discretization [s]
    )
    """
    Assumes an average profile when expanding a time series.
    Thus, the profile is changing in each step of the given time interval of the model.
    """
    N_ave = Int64(dt_disc/dt_data)
    val = nom_val.*profile
    val_ave = reduce(vcat, mean(reshape(val, N_ave, :), dims=1))
    return val_ave
end


function create_profile(
    uniform::Bool, # Uniform (true) or average (false) profile
    nom_val::Number, # Nominal value [p.u.]
    profile::Vector{<:Number}, # Nominalized profile
    dt_data::Int, # Resolution of input data [s]
    dt_disc::Int # Time interval discretization [s]
    )
    if uniform
        type, val = "uniform", uniform_profile(nom_val, profile, dt_data, dt_disc)
    else
        type, val = "average", average_profile(nom_val, profile, dt_data, dt_disc)
    end
    return type, val
end



############### Power network structs ###############
## Components EL network
struct Bus
    id::Int # id
    slack::Int # 1:slack bus, 0:otherwise
    function Bus(
        id::Int, # id
        slack::Int # 1:slack bus, 0:otherwise
        )
        return new(id, slack)
    end
end


abstract type DispatchableGenerator end


struct NonGasFiredGenerator <: DispatchableGenerator
    id::Int # id
    bus::Int # Power bus id
    P_min::Number # Minimum power [p.u.]
    P_max::Number # Maximum power [p.u.]
    P_down::Number # Maximmum ramp down rate [p.u.]
    P_up::Number # Maximmum ramp up rate [p.u.]
    C_l::Number # Linear cost coefficient [p.u.]
    C_q::Number # Quadratic cost coefficient [p.u.]
    function NonGasFiredGenerator(
        id::Int, # id
        bus::Int, # Power bus id
        P_min::Number, # Minimum power [MW]
        P_max::Number, # Maximum power [MW]
        P_down::Number, # Minimum ramp down rate [MW/h]
        P_up::Number, # Maximum ramp up rate [MW/h]
        C_l::Number, # Linear cost coefficient [$/MWh]
        C_q::Number, # Quadratic cost coefficient [$/MW^2h]
        S_base::Number, # Power base [MVA]
        Cost_base::Number # Cost base [$]
    )
        P_min, P_max, P_down, P_up =
            P_min/S_base, P_max/S_base, P_down/S_base, P_up/S_base
        C_l, C_q = C_l*S_base/Cost_base, C_q*S_base^2/Cost_base 

        return new(id, bus, P_min, P_max, P_down, P_up, C_l, C_q)
    end
end


struct GasFiredGenerator <: DispatchableGenerator
    id::Int # id
    bus::Int # Power bus id
    P_min::Number # Minimum power [p.u.]
    P_max::Number # Maximum power [p.u.]
    P_down::Number # Minimum ramp down rate [p.u.]
    P_up::Number # Maximum ramp up rate [p.u.]
    node::Int # Node location in the gas system
    Γ::Number # Power conversion factor [p.u.]
    function GasFiredGenerator(
        id::Int, # id
        bus::Int, # Power bus id
        P_min::Number, # Minimum power [MW]
        P_max::Number, # Maximum power [MW]
        P_down::Number, # Minimum ramp down rate [MW/h]
        P_up::Number, # Maximum ramp up rate [MW/h]
        node::Int, # Node location in the gas system
        Γ::Number, # Power conversion factor [kg/(sMW)]
        S_base::Number, # Power base [MVA]
        M_base::Number # Mass flow base [kg/s]
    )

        P_min, P_max, P_down, P_up =
            P_min/S_base, P_max/S_base, P_down/S_base, P_up/S_base
        Γ = Γ/M_base*S_base

        return new(id, bus, P_min, P_max, P_down, P_up, node, Γ)
    end
end


#TODO: Do we need profile type?
struct WindGenerator
    id::Int # id
    bus::Int # Power bus id
    P_max::Number # Maximum power [p.u.]
    type_profile::AbstractString # Type of wind load profile
    profile_type::AbstractString # Type of wind profile
    wind_profile::Vector{Float64} # Normalized daily wind profile
    P_w::Vector{Float64} # Wind power forecast in each timestep [p.u.]
    function WindGenerator(
        id::Int, # id
        bus::Int, # Power bus id
        P_max::Number, # Maximum power [MW]
        profile_type::AbstractString, # Type of wind profile
        wind_profile::Vector{Float64}, # Normalized daily wind profile
        dt::Int64, # Seconds per timestep in model
        dt_data::Number, # Seconds per timestep in data
        uniform::Bool, # Uniform (true) or average (false) profile
        S_base::Number # Power base [MVA]
        )

        P_max = P_max/S_base
        type, P_w = create_profile(uniform, P_max, wind_profile, dt_data, dt)

        return new(id, bus, P_max, type, profile_type, wind_profile, P_w)
    end
end


struct PowerDemand
    id::Int # id
    bus::Int # Power bus id
    load_nom::Number # Nominal electrical load [p.u.]
    type_profile::AbstractString # Type of power load profile
    load_profile::Vector{Float64} # Electrical load profile
    D_el::Vector{Float64} # Electrical load at each timestep
    function PowerDemand(
        id::Int, # id
        bus::Int, # Power bus id
        load_nom::Number, # Nominal electrical load [MW]
        load_profile::Vector{Float64}, # Electrical load profile
        dt::Int64, # Time interval discretization [s]
        dt_eload_s::Number, # Resolution of input data [s]
        uniform::Bool, # Uniform (true) or average (false) profile
        S_base::Number # Power base [MVA]
        )

        load_nom = load_nom/S_base
        type, D_el = create_profile(uniform, load_nom, load_profile, dt_eload_s, dt)

        return new(id, bus, load_nom, type, load_profile, D_el)
    end
end


struct Line
    id::Int # id
    start::Int # Powerline start
    stop::Int # Powerline ends
    X::Float64 # Line reactance in [p.u.]
    F_max::Number # Maximum capacity [p.u.]
    ptdf::Vector{Float64} # PTDFs [dimensionless]

    function Line(
        id::Int, # id
        start::Int, # Powerline start
        stop::Int, # Powerline ends
        X::Float64, # Line reactance in [p.u.]
        F_max::Number, # Maximum capacity [MW]
        ptdf::Vector{Float64}, # PTDFs [dimensionless]
        S_base::Number # Power base [MVA]
    )

        F_max = F_max/S_base
        return new(id, start, stop, X, F_max, ptdf)
    end
end
