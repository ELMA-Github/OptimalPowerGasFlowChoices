cd(dirname(@__FILE__))
using Pkg
Pkg.activate("")

using JuMP, Ipopt, Gurobi, CSV, DataFrames, LinearAlgebra, Dates, XLSX, Statistics
#using Plots
include("components.jl")
include("bases_choice.jl")
include("discretization.jl")
include("input_data.jl")
include("algorithms.jl")
include("energy_system.jl")
include("PTDF_matrix.jl")
include("variables.jl")
include("constraints.jl")
include("objective.jl")
include("model_wrapper.jl")
include("save_results.jl")
#include("spatial_branch_and_bound.jl")
# cd(dirname(@__FILE__))

config_dict = Dict{Symbol,Any}(
    ########## System must be either 'integrated_power_and_gas' or 'gas_only'
    :system => "integrated_power_and_gas",
    # Modeltype: choose one of
    # 'NCNLP', 'SCP', 'MINLP', 'MILP', 'MISOCP', 'MISOCP_McCormick', 'MISOCP_weymouth', 'MISOCP_sbnb', 'PWL', 'CELP', trade'
    :model_type => "NCNLP",
    ########## Choose objective ##########
    # One of "linear" or "quadratic"
    :objective_type => "quadratic",
    ########## Choose warmstart ##########
    # If you don't want to warmstart, pass an empty Dict
    # :warmstart => Dict{Symbol,Any}(:method => "MILP"),#
    :warmstart => Dict{Symbol,Any}(),
    ########## Choose case study ##########
    # Gas_tree: 4-node, 3 pipes, 1 source, uniform load #(provided every hour in the input data)
    # Gas_tree_2sources: 4-node, 3 pipes, 2 surce, uniform load #(provided every hour in the input data)
    # Gas_tree_2sources_3min: 4-node, 3 pipes, 2 source, # average load (provided every 3 min in the input data)
    # :case_study => "Gas_line", "Case_Study_A", "GasLib-40_IEEE24"
    :case_study => "GasLib-40_IEEE24",
    ########## Time discretization ##########
    :timehorizon => 24, # Time horizon (hours)
    :dt => 900, # Time interval for discretization (seconds)
    ########## Space discretization ##########
    # 0: No space discretization --> discretization segments equal pipeline length
    # 1: Space discretization--> user defines each segments' length
    :space_disc => 1,
    :segment => 50000, #only relevant if space_disc=1
    ########## Per unit conversion ##########
    :pu => true,
    ########## Other Parameters ##########
    :print_model => 0, # Control for model printing
    :C_cur_load_gas => 36000, # Value of lost load ($/(kgh/s))
    :C_cur_load_power => 1000, # Value of lost load ($/MWh)
    # :M => 10000 # Parameter for Big-M in MISCOP
)

if config_dict[:system] == "gas_only"
    pop!(config_dict, :C_cur_load_power, 0)
end

sys = config_dict[:system]
config_dict_algorithm = Dict{Symbol,Any}()
########## Select additional model parameters for various model formulations
if config_dict[:model_type] == "PWL"
    ########## Parameters PWL ##########
    # Number of set points for mass flow and average pressure {Int}
    config_dict_algorithm[:no_m_set] = 3 # >= 3
    config_dict_algorithm[:no_π_set] = 3 # >= 3
    config_dict_algorithm[:method] = "ldcc" # One of ["ldcc"]
elseif config_dict[:model_type] == "SCP"
    ########## Parameters SCP ##########
    config_dict_algorithm[:order] = 1 # Taylor-series expansion order, one of {1,2}
    config_dict_algorithm[:k_max] = 100 # Maximum number of iterations
    config_dict_algorithm[:δ_init] = 10e-3
    config_dict_algorithm[:δ_update] = 2
    config_dict_algorithm[:δ_max] = 10e3
    config_dict_algorithm[:ϕ_max] = 10e-9
elseif config_dict[:model_type] == "MISOCP_sbnb"
    config_dict_algorithm[:cut_slack] = 0
    config_dict_algorithm[:max_relaxation_gap] = 0.1
    config_dict_algorithm[:max_optimal_physics_gap] = 0.1
    config_dict_algorithm[:abort] = false
    # config_dict_algorithm[:M] = 5000
elseif config_dict[:model_type] in ["MILP", "MISOCP"]
    config_dict_algorithm[:linear_overestimator] = false
end

if config_dict[:model_type] in ["SCP", "NCNLP", "MINLP", "MISOCP","MISOCP_McCormick", "MISOCP_sbnb", "CELP", "PWL"]
    ########## Parameters for PDE-based models ##########
    # (1: steady state, 2: quasi-dynamic, 3: transient) 
    config_dict_algorithm[:PDE_type] = 3
end

# if config_dict[:model_type] in ["MISOCP_weymouth"]
#     ########## Parameters for PDE-based models ##########
#     # Attention! This is set by definition! 
#     config_dict_algorithm[:PDE_type] = 2
# end

# # At least the solver name MUST be set
# config_dict_solver = Dict{Symbol,Any}(
#     :solver_name => "Gurobi",
#     # :max_cpu_time => 100.0,
#     # :NonConvex => 2,
#     :NumericFocus => 3,
#     :TimeLimit => 1000,
#     # :BarConvTol => 1e-4,
#     # :FeasibilityTol => 1e-9
#     # :MIPFocus => 1
#     )
        

# boundary_conditions = get_boundary_conditions(
#     deepcopy(config_dict), 3)

# ES = intialize_energy_system(config_dict)
# intialize_data!(ES);
# ES.boundary_conditions = Dict(
#     :m => boundary_conditions[:m]/ES.bases_dict[:M_base],
#     :π_avg => boundary_conditions[:π_avg]/ES.bases_dict[:Π_base],
# )

# ## start = Dates.now()
# run_model!(ES, config_dict_algorithm)#; config_dict_solver)
# ## stop = Dates.now()
# ## println("Solution time: $((stop-start)/Millisecond(1000)) seconds.")



config_dict[:dt] = 900
boundary_conditions = get_boundary_conditions(
    deepcopy(config_dict), 3)

model_list = ["NCNLP", "SCP", "CELP", "MISOCP", "MILP"]
LO = [true,false]
PDE_types = [3]
dt_list = [900]#, 300]#[3600, 900]#, 900]
results = Dict{String, EnergySystem}()
for dt in dt_list

    for PDE_type in PDE_types

        config_dict[:dt] = dt
        if dt in [300,900]
            config_dict[:space_disc] = 0
            config_dict[:segment] = 50000
        else
            config_dict[:space_disc] = 0
        end
        config_dict[:system] = "integrated_power_and_gas" 


        ES = intialize_energy_system(config_dict)
        intialize_data!(ES);
        ES.boundary_conditions = Dict(
            :m => boundary_conditions[:m]/ES.bases_dict[:M_base],
            :π_avg => boundary_conditions[:π_avg]/ES.bases_dict[:Π_base],
        )

        for model in model_list

            ES.model_type = model
            config_dict[:system] = "integrated_power_and_gas"

            if (model == "MISOCP") # && (dt == 900)
                ES.warmstart = Dict{Symbol,Any}(:method => "MILP")
            elseif model == "SCP"
                ES.warmstart = Dict{Symbol,Any}(:method => "CELP")
            else
                ES.warmstart = Dict{Symbol,Any}()
            end

            config_dict_algorithm = Dict{Symbol,Any}()
            ########## Select additional model parameters for various model formulations
            if model == "PWL"
                ########## Parameters PWL ##########
                # Number of set points for mass flow and average pressure {Int}
                config_dict_algorithm[:no_m_set] = 3 # >= 3
                config_dict_algorithm[:no_π_set] = 3 # >= 3
                config_dict_algorithm[:method] = "ldcc" # One of ["ldcc"]
            elseif model == "SCP"
                ########## Parameters SCP ##########
                config_dict_algorithm[:order] = 1 # Taylor-series expansion order, one of {1,2}
                config_dict_algorithm[:k_max] = 100 # Maximum number of iterations
                config_dict_algorithm[:δ_init] = 10e-3
                config_dict_algorithm[:δ_update] = 2
                config_dict_algorithm[:δ_max] = 10e3
                config_dict_algorithm[:ϕ_max] = 10e-12
            end

            if model in ["SCP", "NCNLP", "MINLP", "MISOCP","MISOCP_McCormick", "MISOCP_sbnb", "CELP", "PWL"]
                ########## Parameters for PDE-based models ##########
                config_dict_algorithm[:PDE_type] = PDE_type
            end

            if model in ["MISOCP", "MILP"]
                for lo in LO
                    config_dict_algorithm[:linear_overestimator] = lo
                    run_model!(ES, config_dict_algorithm)
                    lo_string = lo == false ? "_"*string(false) : ""
                    results[model*lo_string*"_"*string(PDE_type)*"_"*string(dt)] = deepcopy(ES)
                end
            else
                run_model!(ES, config_dict_algorithm)
                results[model*"_"*string(PDE_type)*"_"*string(dt)] = deepcopy(ES)
            end

        end
    end
end


# # save_results_to_excel(ES, config_dict)
