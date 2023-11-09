"""
Created by Enrica Raheli and Yannick Werner
Licensed under xxxx
Language: Julia, version number
----------------------------------
Content of this file:
Function definitions for space discretization of a gas network.
"""

function segment_discretization_equal(L, segment_length)
    """
    Divide pipe in equal segments of given length divided by the number of segments created.
    """
    N_dx = ceil(L/segment_length) # number of segments of equal length
    seg_length = L/N_dx # length of each segment
    last_length = seg_length
    return N_dx, seg_length, last_length
end

function segment_discretization_maximum(L, segment_length)
    """
    Divide pipe in equal segments of given length and one potentially smaller segment.
    """
    N_seg = div(L, segment_length) # number of segments of constant length
    if rem(L, segment_length) == 0 # remainder euqual to zero -> length is a multiple of discretization segment
        N_dx = N_seg # Number of segments in the pipe = number of segments of constant lenght
        last_length = segment_length # length of the last subpipe
    else
        N_dx = N_seg+1 # Number of segments in the pipe = number of segments of constant lenght + one smaller segment
        last_length = rem(pipe.L, segment_length); # length of the last subpipe
    end
    seg_length = segment_length
    return N_dx, seg_length, last_length
end


function space_discretization(nodes, pipes, segment, bases_dict::Dict{Symbol,Number}, type="equal")
    """
    Divide pipe in equal segments of given length and one potentially smaller segment.
    """
    @info "Start space discretization."
    node_id_count, pipe_id_count = maximum(keys(nodes))+1, maximum(keys(pipes))+1
    sub_nodes = nodes
    sub_pipelines = Dict{Int64, Pipeline}()
    segment_length = segment/bases_dict[:D_base]

    if type == "equal"
        space_discretization_method = segment_discretization_equal
    elseif type == "maximum"
        space_discretization_method = space_discretization_maximum
    else
        throw(ArgumentError(
            "Space discretization type must be one of 'equal' or 'maximum'."))
    end

    for (id, pipe) in pipes

        N_dx, seg_length, last_length = space_discretization_method(pipe.L, segment_length)

        # Calculate number of auxiliary subnodes between original nodes
        N_nodes_between = Int64(N_dx - 1)
        
        x_coord_nodes = LinRange(nodes[pipe.start].x, nodes[pipe.stop].x, N_nodes_between+2)
        y_coord_nodes = LinRange(nodes[pipe.start].y, nodes[pipe.stop].y, N_nodes_between+2)
    
        if N_nodes_between == 0
            sub_pipelines[id] = pipe
        else
            pipe_sub_nodes = [] # subnodes of discretized pipe
            for n in 1:N_nodes_between
                ### Creation of (sub)nodes
                # Note that this covers the case of parallel pipes to not create multiple
                # subnodes at the same location.
                node_id = node_id_count #parse(Int64, string((pipe.start),(pipe.stop),n))
                sub_nodes[node_id] = create_subnode(nodes, pipe, node_id, x_coord_nodes[n+1], y_coord_nodes[n+1])
        
                push!(pipe_sub_nodes, node_id)
                ### Creation of (sub)pipelines
                if n == 1
                    # parse(Int64, string(id,(pipe.start),(pipe.stop),n))
                    sub_pipelines[id] = create_subpipe(
                        pipe, id, id, seg_length, sub_nodes, pipe.start, pipe_sub_nodes[n], bases_dict)
                else
                    sub_pipelines[pipe_id_count] = create_subpipe(
                        pipe, pipe_id_count, id, seg_length, sub_nodes, pipe_sub_nodes[n-1], pipe_sub_nodes[n], bases_dict)
                    pipe_id_count += 1
                end
                
                if n == N_nodes_between
                    sub_pipelines[pipe_id_count] = create_subpipe(
                        pipe, pipe_id_count, id, last_length, sub_nodes, pipe_sub_nodes[n], pipe.stop, bases_dict)
                    pipe_id_count += 1
                end

                node_id_count += 1
            end
        end
    end
    @info "End space discretization. $(length(sub_pipelines)-length(pipes)) new pipelines created."
    return sub_nodes, sub_pipelines
end


function create_subnode(nodes, pipe, node_id, x, y)
    """
    Creates auxiliary subnodes between two original nodes 
    if the pipeline between them is discretized.
    """

    # The new subnodes don't have explicit pressure limits.
    # Hence, the minimum and maximum are used as dummy here.
    pmax_sub = max(nodes[pipe.start].Π_max,nodes[pipe.stop].Π_max)
    pmin_sub = min(nodes[pipe.start].Π_min,nodes[pipe.stop].Π_min)

    return NonSlackNode(node_id, pmax_sub, pmin_sub, x, y, 1)
end


function create_subpipe(
    pipe,
    subpipeline_id,
    parentpipe_id,
    pipe_length,
    sub_nodes,
    start_node_id::Int,
    end_node_id::Int,
    bases_dict::Dict{Symbol,Number})
    """
    Creates auxiliary subpipeline between two (sub-)nodes.
    """

    D_base, C_base = bases_dict[:D_base], bases_dict[:C_base]

    return Pipeline(subpipeline_id, start_node_id, end_node_id, parentpipe_id,
        pipe.F, pipe.D*D_base, pipe_length*D_base, sub_nodes, D_base, C_base)

end