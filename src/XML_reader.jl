using EzXML, XMLDict, DataFrames, OrderedCollections, CSV
using Distributions

cd(dirname(@__FILE__))

function convert_dict_types(dict::OrderedDict{Any, Any})
    return Dict(Symbol(replace(String(key), "-" => "")) => val for (key,val) in dict)
end


function get_unit(val::String)
    return ""
end


function get_unit(val::OrderedDict{Any, Any})
    return get(val, :unit, "")
end


function get_unit_df(dict::Union{Dict{Symbol, Any}, OrderedDict{Any, Any}})
    unit_dict = Dict(key => get_unit(val) for (key,val) in dict)
    unit_df = DataFrame(unit_dict)
    return unit_df
end



function get_value(val::String)
    return val
end


function get_value(val::OrderedDict{Any, Any})
    return get(val, :value, "")
end


function get_value_df(dict::Dict{Symbol, Any})
    val_dict = Dict(key => get_value(val) for (key,val) in dict)
    val_df = DataFrame(val_dict)
    return val_df
end


function create_data_df(data_dict::OrderedDict{Any, Any}, node_type::AbstractString)

    if isa(data_dict[node_type], Vector{Any})
        data_dict = convert_dict_types.(data_dict[node_type]) 
        unit_df = get_unit_df(data_dict[1])
        val_df = reduce(vcat, get_value_df.(data_dict))
    else
        data_dict = convert_dict_types(data_dict[node_type]) 
        unit_df = get_unit_df(data_dict)
        val_df = get_value_df(data_dict)
    end
    return unit_df, val_df
end


function read_xml(path::AbstractString)

    doc = readxml(path)
    xml_data = xml_dict(doc)
    data = xml_data["network"]

    model_data = Dict()
    model_units = Dict()

    data_nodes = data["nodes"]
    for node_type in keys(data_nodes)
        unit_df, val_df = create_data_df(data_nodes, node_type)
        model_data[Symbol(node_type)] = val_df
        model_units[Symbol(node_type)] = unit_df
    end

    data_connections = data["connections"]
    for node_type in keys(data_connections)
        unit_df, val_df = create_data_df(data_connections, node_type)
        model_data[Symbol(node_type)] = val_df
        model_units[Symbol(node_type)] = unit_df
    end

    return model_data, model_units
end


function friction_factor_colebrook_white(
    e::Number, # roughness length
    D::Number, # inner pipe diameter
    )
    """
    E. Shashi Meno (2005). Gas pipeline hydraulics, pp.47--50.
    Assuming large Reynold numbers, turbulent flow, fully rough.
    """
    return (-2*log10(e/(3.7*D)))^(-2)
end


function transform_df_string_to_float!(df)
    for col in names(df)
        try transform!(df, col => ByRow(x -> parse(Float64, x)) => col)
        catch ArgumentError
        end
    end
    return nothing
end


function F1(z::Number,b::Vector{<:Number})
    return b'*[1,z,z^2]
end


function F2(x::Number,y::Number,A::Matrix{<:Number})
    return [1,x,x^2]'*A*[1,y,y^2]
end

function load_mapping(path, n, type)
    path = path*"GasLib-$(string(n))-v1-20211130_raw///mapping_$type.csv"
    file = CSV.read(path, DataFrame)
    try select!(file, Not([:from_id_type, :to_id_type])) catch e end
    # deletecols!(file, :from_id_type)
    # deletecols!(file, :to_id_type)
    return file
end

rename_dict_nodes = Dict(
    "id" => "id_type",
    "pressureMax" => "Pmax_MPa",
    "pressureMin" => "Pmin_MPa",
    "x" => "x",
    "y" => "y"
)

rename_dict_pipes = Dict(
    "length" => "Length_m",
    "diameter" => "Diameter_m",
    "id" => "id_type",
    "friction" => "friction",
    "from" => "from_id_type",
    "to" => "to_id_type"
)

rename_dict_compressors = Dict(
    "pressureInMin" => "Pmin_in_MPa",
    "pressureOutMax" => "Pmax_out_MPa",
    "id" => "id_type",
    "from" => "from_id_type",
    "to" => "to_id_type",
    "fuelGasVertex" => "fuel_gas_node"
)

no_nodes_case_study = [40]#[11, 135]
node_types = [:source, :innode, :sink]
data_sets = Dict()
unit_sets = Dict()
data_final = Dict()
path_data = "..//inputs//gas_only//"

for n in no_nodes_case_study
    path = path_data*"GasLib-$(string(n))-v1-20211130_raw//GasLib-$(string(n))-v1-20211130.net"
    data_sets[n], unit_sets[n] = read_xml(path)
end

# normDensity, molarMass, calorificValue = 
#     data_sets[11][:source][1, ["normDensity", "molarMass", "calorificValue"]] # kg_per_m_cube, kg_per_kmol, MJ_per_m_cube
# normDensity, molarMass, calorificValue = parse.(Float64, [normDensity, molarMass, calorificValue])

const SPEED_OF_SOUND = 350
# function weymouth_coeff(D, F, L)
# K = sqrt(D*A^2/(F*C^2*L))

for n in no_nodes_case_study
    data_final[n] = Dict()

    ############### Nodes ###############
    for n_type in node_types
        select!(data_sets[n][n_type], keys(rename_dict_nodes)...)
        rename!(data_sets[n][n_type], rename_dict_nodes)
        # delete!(data_sets[n][n_type], 1)
        transform_df_string_to_float!(data_sets[n][n_type])
    end

    global nodes_df = vcat([data_sets[n][n_type] for n_type in node_types]...)
    select!(nodes_df, Not([:x,:y]))

    # nodes_df[!,"Node_No"] = 1:nrow(nodes_df)
    leftjoin!(nodes_df, load_mapping(path_data, n, "nodes"), 
        on="id_type")

    # data_final[n]["nodes"] = nodes_df
    for n_type in node_types
        leftjoin!(data_sets[n][n_type], nodes_df[!,["id_type", "Node_No"]], on = "id_type")
        rename!(data_sets[n][n_type], "Node_No" => "Node")
    end

    if n == 40
        slack_node = ["source_2", "source_3"]
    elseif n == 11
        slack_node = ["entry01"]
    end

    slack_node_bin = [row.id_type in slack_node for row in eachrow(nodes_df)]

    nodes_df[!, "Node_Type"] .= 0
    nodes_df[slack_node_bin, "Node_Type"] .= 1
    nodes_df[!, "Pslack_MPa"] .= NaN
    nodes_df[slack_node_bin, "Pslack_MPa"] .= 
        nodes_df[slack_node_bin, "Pmax_MPa"]/1.5

    sort!(nodes_df, :Node_No)

    sort!(data_sets[n][:source], :Node)
    data_sets[n][:source][!,"Supply_No"] = 1:nrow(data_sets[n][:source])
    # data_sets[n][:source][!,"Smax_kg_s"] .= 725*0.785/3.6 #2000.0
    data_sets[n][:source][!,"Smin_kg_s"] .= 0.0
    # d_C1, d_C2 = Normal(0.15,0.015), Normal(0.01,0.001)
    # data_sets[n][:source][!,"C1_per_kgh"] = rand(d_C1, nrow(data_sets[n][:source]))
    # data_sets[n][:source][!,"C2_per_kgh2"] = rand(d_C2, nrow(data_sets[n][:source]))

    load_data = CSV.read(
        path_data*"GasLib-$(string(n))-v1-20211130_raw///sources_attributes.csv", DataFrame)
    leftjoin!(data_sets[n][:source], load_data, on="id_type")

    sort!(data_sets[n][:sink], :Node)
    data_sets[n][:sink][!,"Load_No"] = 1:nrow(data_sets[n][:sink])
    # data_sets[n][:sink][!,"Load_kg_s"] .= 75*0.785/3.6 #7.0
    # data_sets[n][:sink][!,"Profile"] .= "Gas_profileA"
    load_data = CSV.read(
        path_data*"GasLib-$(string(n))-v1-20211130_raw///sinks_attributes.csv", DataFrame)
    leftjoin!(data_sets[n][:sink], load_data, on="id_type")

    ############### Pipes ###############
    global pipe_df = deepcopy(data_sets[n][:pipe])

    if n == 11
        push!(pipe_df, pipe_df[1,:])
        pipe_df[9,[:from,:to,:id]] = ["N01", "N03", "pipe09_N01_N03"]
    end


    transform_df_string_to_float!(pipe_df)

    pipe_df[!,"length"] *= 1000 # convert to m
    pipe_df[!,"diameter"] /= 1000 # convert to m
    pipe_df[!,"roughness"] /= 1000 # convert to m
    pipe_df[!, "friction"] = friction_factor_colebrook_white.(
        pipe_df[!, "roughness"], pipe_df[!, "diameter"])

    select!(pipe_df, keys(rename_dict_pipes)...)
    rename!(pipe_df, rename_dict_pipes)

    leftjoin!(pipe_df, nodes_df[!,["id_type", "Node_No"]], 
        on="from_id_type"=>"id_type")
    rename!(pipe_df, "Node_No" => "From_Node")
    leftjoin!(pipe_df, nodes_df[!,["id_type", "Node_No"]], 
        on="to_id_type"=>"id_type")
    rename!(pipe_df, "Node_No" => "To_Node")

    leftjoin!(pipe_df, load_mapping(path_data, n, "pipes"), on="id_type")
    # pipe_df[!,"Pipe_No"] = 1:nrow(pipe_df)
    for row in eachrow(pipe_df)
        if ~ismissing(row.switch_dir) && (row.switch_dir == 1)
            from = row.From_Node
            row.From_Node = row.To_Node
            row.To_Node = from
        end
    end
    sort!(pipe_df, :Pipe_No)

    ############### Compressors ###############
    global compressor_df = data_sets[n][:compressorStation]
    transform_df_string_to_float!(compressor_df)
    select!(compressor_df, keys(rename_dict_compressors)...)
    rename!(compressor_df, rename_dict_compressors)
    leftjoin!(compressor_df, nodes_df[!,["id_type", "Node_No"]], 
        on="from_id_type"=>"id_type")
    rename!(compressor_df, "Node_No" => "From_Node")
    leftjoin!(compressor_df, nodes_df[!,["id_type", "Node_No"]], 
        on="to_id_type"=>"id_type")
    rename!(compressor_df, "Node_No" => "To_Node")

    leftjoin!(compressor_df, nodes_df[!,["id_type", "Node_No"]], on = [:fuel_gas_node => :id_type])
    select!(compressor_df, Not(:fuel_gas_node))
    rename!(compressor_df, :Node_No => :fuel_gas_node)
    compressor_df[!, :fuel_gas_consumption] .= 0.005

    # compressor_df[!,"Compressor_No"] .= 1:nrow(compressor_df)
    leftjoin!(compressor_df, load_mapping(path_data, n, "compressors"), on="id_type")
    compressor_df[!,"Compression_cost"] .= 2.0
    compressor_df[!,"CR_Max"] .= 1.5
    compressor_df[!,"CR_Min"] .= 1.0

    # Limit pressure on nodes before a compressor
    for row in eachrow(compressor_df)
        nodes_df[nodes_df[!,"id_type"] .== row.from_id_type, "Pmin_MPa"] .= row.Pmin_in_MPa
    end

    # Limit logically pressure at nodes before and after compressor
    # if n == 40
    #     above_min = ["sink_13", "sink_14", "sink_10", "sink_11", "innode_3", "innode_2", "sink_19", "sink_27", "innode_5", "innode_4", "sink_23", "source_1"]
    #     not_below_max = ["source_2", "sink_3", "sink_23", "source_1", "source_3"]
    #     flt_min =  [row.id_type in above_min for row in eachrow(nodes_df)]
    #     nodes_df[flt_min, "Pmin_MPa"] .= 31.0132
    #     flt_max =  [~(row.id_type in not_below_max) for row in eachrow(nodes_df)]
    #     nodes_df[flt_max, "Pmax_MPa"] .= 71.0132
    # end

    # Convert units
    nodes_df[!, ["Pmax_MPa", "Pmin_MPa", "Pslack_MPa"]] ./= 10 # bar to MPa
    nodes_df[!, "Pmin_MPa"] .= 3.101325

    ###### Remove nodes and pipes ######
    if n == 40
        filter!(row -> row.id_type != "innode_3", nodes_df)
        filter!(row -> ~(row.id_type in ["pipe_33", "pipe_39"]), pipe_df)
    end

    # Save files
    path_save = path_data*"GasLib-$(string(n))"
    if isdir(path_save)
        nothing
    else
        mkdir(path_save)
    end

    # data_final[n]

    col_order_pipes = ["Pipe_No", "From_Node", "To_Node", "Length_m", "Diameter_m", "friction"]
    col_order_nodes = ["Node_No", "Pmin_MPa", "Pmax_MPa", "Pslack_MPa", "Node_Type", "x", "y"]
    col_order_supply = ["Supply_No", "Node", "Smax_kg_s", "Smin_kg_s", "C1_per_kgh", "C2_per_kgh2"]
    col_order_load = ["Load_No", "Node", "Load_kg_s", "Profile"]
    col_order_compressor = ["Compressor_No","From_Node", "To_Node", "fuel_gas_node", "fuel_gas_consumption", 
        "CR_Max", "CR_Min", "Compression_cost"]
    CSV.write(joinpath(path_save,"gas_pipes.csv"), pipe_df[!, col_order_pipes])
    CSV.write(joinpath(path_save,"gas_supply.csv"), data_sets[n][:source][!, col_order_supply])
    CSV.write(joinpath(path_save,"gas_load.csv"), data_sets[n][:sink][!, col_order_load])
    CSV.write(joinpath(path_save,"gas_nodes.csv"), nodes_df[!, col_order_nodes])
    CSV.write(joinpath(path_save,"gas_compressors.csv"), compressor_df[!, col_order_compressor])

end 



const T_amb = 293.15

n = 40
path_compressor_data = path_data*"GasLib-$(string(n))-v1-20211130_raw//GasLib-$(string(n))-v1-20211130.cs"
function read_xml_compressor(path::AbstractString)

    doc = readxml(path)
    xml_data = xml_dict(doc)
    data_compressors = xml_data["compressorStations"]["compressorStation"]

    return data_compressors
end

# #     model_data = Dict()
# #     model_units = Dict()

data_compressors = read_xml_compressor(path_compressor_data)

    # for compStation in data_compressors
    #     nomSpeed_per_min = parse(Float64, 
    #         compStation["configurations"]["configuration"]["stage"]["compressor"][:nominalSpeed])
    #     nomSpeed_per_s = nomSpeed_per_min/60

    #     ############### Compressor ###############



    #     ############### Drive ###############
    #     drive_dict = compStation["drives"]["gasTurbine"]
    #     power_fun_coeff_mat = Matrix(#transpose(
    #         parse.(Float64, reshape(
    #         [drive_dict["power_fun_coeff_$a"][:value] for a in 1:9], 3, 3)))#)#
    #     power_fun_coeff_vec = parse.(Float64, [drive_dict["power_fun_coeff_$a"][:value] for a in 1:3])
    #     P_max_drive = F1(nomSpeed_per_s, power_fun_coeff_vec)
    #     P_max_drive = F2(nomSpeed_per_min, T_amb, power_fun_coeff_mat)




#     if isa(data, Vector{Any})
#         data_dict = convert_dict_types.(data) 
#         unit_df = get_unit_df(data_dict[1])
#         val_df = reduce(vcat, get_value_df.(data_dict))
#     else
#         data_dict = convert_dict_types(data) 
#         unit_df = get_unit_df(data_dict)
#         val_df = get_value_df(data_dict)
#     end

#     # unit_df, val_df = create_data_df(data_nodes, node_type)
#     model_data[Symbol(node_type)] = val_df
#     model_units[Symbol(node_type)] = unit_df
#     # end

#     return data
# end

# cs_data = read_xml_compressor(path_data*"GasLib-$(string(n))-v1-20211130//GasLib-$(string(n))-v1-20211130.cs")