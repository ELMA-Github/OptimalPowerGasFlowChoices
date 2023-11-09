function SBNB(
    ES::EnergySystem,
    config_dict_solver::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    config_dict_algorithm::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
    # cut_slack::Number=0,
    # max_relaxation_gap::Number=0.01,
    # PDE_type::Int=3,
    # M::Real=1000
    )
    """
    MISOCP with spatial branch and bound
    For the callback to work, three settings are necessary:
    - use of direct_model
        -> such that the underlying model can be manipulated during optimization
    - CPXPARAM_Threads = 1
        -> to avoid communication failures between threads when using custom callback
    - CPX_PARAM_PREIND = 0
        -> to ensure that the mapping between JuMP and CPLEX models is not changed during presolve
    """

    if ~isempty(ES.warmstart)
        warmstart!(ES, ES.warmstart[:method])
    end

    cut_slack = pop!(config_dict_algorithm, :cut_slack, 0)
    max_relaxation_gap = pop!(config_dict_algorithm, :max_relaxation_gap, 0.01)
    max_optimal_physics_gap = pop!(config_dict_algorithm, :max_optimal_physics_gap, 0.1)
    
    # max_relaxation_gap::Number=0.01,
    # PDE_type::Int=3,
    # M::Real=1000

    global model = direct_model(CPLEX.Optimizer())
    build_model(model, ES, config_dict_algorithm; model_type="MISOCP_sbnb")

    if ~isempty(ES.warmstart)
        set_warmstart!(model, "MISOCP_sbnb", ES)
    end
    
    ES.model = model

    pipes, P, T = ES.pipes, ES.P, ES.T

    @variables model begin
        pipes[p].Π_avg_min <= π_avg_b[p in P, t in T] <= pipes[p].Π_avg_max
        0 <= γ⁺_b[p in P, t in T] <= pipes[p].Γ_max
        0 <= γ⁻_b[p in P, t in T] <= pipes[p].Γ_max
    end

    # set_solver(model, model_type, config_dict_solver)
    @info "Setting CPLEX attribute Threads = 1 and PREIND (presolve) = 0
    to use custom branching rules."
    set_optimizer_attribute(model, "CPXPARAM_Threads", 1) 
    set_optimizer_attribute(model, "CPX_PARAM_PREIND", 0)
    set_optimizer_attribute(model, "CPX_PARAM_TILIM", 100)
    # set_optimizer_attribute(model, "CPX_PARAM_BARQCPEPCOMP", 10^-4)
    set_optimizer_attribute(model, "CPX_PARAM_BARQCPEPCOMP", 10^-3)
    # set_optimizer_attribute(model, "CPX_PARAM_EPINT", 10^-3)

     
    # set_optimizer_attribute(model, "CPX_PARAM_EPRHS", 1e-08)
    # set_optimizer_attribute(model, "CPX_PARAM_BAREPCOMP", 1e-08)
    

    set_optimizer_attribute(model, "CPX_PARAM_MIQCPSTRAT", 1)

    Γ_max = [pipes[p].Γ_max for p in P, t in T] 
    slope_max_point = [pipes[p].Γ_max/pipes[p].Π_avg_max
        for p in P, t in T] 

        # CPX_PARAM_MIQCPSTRAT, 1

    # set_optimizer_attribute(model, "CPX_PARAM_CLIQUES", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_FPHEUR", -1)
    # set_optimizer_attribute(model, "CPXPARAM_MIP_Cuts_LiftProj", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_LANDPCUTS", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_COVERS", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_DISJCUTS", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_FLOWCOVERS", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_FLOWPATHS", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_FRACCUTS", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_GUBCOVERS", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_ZEROHALFCUTS", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_MIRCUTS", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_IMPLBD", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_PROBE", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_HEURFREQ", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_RINSHEUR", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_CUTPASS", -1)
    # set_optimizer_attribute(model, "CPX_PARAM_PREIND", 0)
    # set_optimizer_attribute(model, "CPX_PARAM_AGGIND", 0)
    # set_optimizer_attribute(model, "CPX_PARAM_RELAXPREIND", 0)
    # set_optimizer_attribute(model, "CPX_PARAM_PREPASS", 0)
    # set_optimizer_attribute(model, "CPX_PARAM_REPEATPRESOLVE", 0)
    # set_optimizer_attribute(model, "CPX_PARAM_BNDSTRENIND", 0)
    # set_optimizer_attribute(model, "CPX_PARAM_COEREDIND", 0)
    # set_optimizer_attribute(model, "CPX_PARAM_CUTSFACTOR", 1)
    
    # set_optimizer_attribute(model, "CPX_PARAM_MIPEMPHASIS", 1)
    # set_optimizer_attribute(model, "CPX_PARAM_HEURFREQ", -1)
    # set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Reformulations", 0)


    if ~isempty(config_dict_solver)
        delete!(config_dict_solver, :solver_name)
        for (key, val) in config_dict_solver
            set_optimizer_attribute(model, String(key), val)
        end
    end

    obj_val_non_convex_candidate = 1e75
    min_max_relaxation_gap = 1e75
    sol = ()
    
    function my_callback_function(cb_data::CPLEX.CallbackContext, context_id::Clong)

        @info context_id

        if context_id != CPX_CALLBACKCONTEXT_CANDIDATE && context_id != CPX_CALLBACKCONTEXT_BRANCHING
            return
        end

        @info "start callback"

        node_id = Ref{Cint}()
        ret = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODEUID, node_id)
        @info "Node id: $(node_id[])"

        x_p, obj_p = Vector{Cdouble}(undef, 1), Ref{Cdouble}() # Define empty containeers 

        if context_id == CPX_CALLBACKCONTEXT_CANDIDATE

            ret = CPXcallbackgetcandidatepoint(cb_data, x_p, 0, 0, obj_p)
            obj_p = obj_p[]
            println("Current candidate point objective: $obj_p.")

            # ispoint = Ref{Cint}()
            # CPXcallbackcandidateispoint(cb_data, ispoint)
            # println("Current candidate is integer feasible: $(ispoint[]).")

            # if (obj_val_non_convex_candidate - obj_p) / obj_val_non_convex_candidate <= max_optimal_physics_gap
            if obj_p > obj_val_non_convex_candidate*1.001
                @info "Reject candidate because of optimality gap: $obj_p <= $obj_val_non_convex_candidate."
                ret = CPXcallbackrejectcandidate(
                    cb_data, Cint(0), Cint(0), Cdouble[], Ref{Cchar}(), Ref{Cint}(), Cint[], Cdouble[])
                return nothing
            # elseif obj_p <= obj_val_non_convex_candidate
            #     @info "Update best bound."
            #     obj_val_non_convex_candidate = obj_p
            end
            
        end

        if context_id == CPX_CALLBACKCONTEXT_BRANCHING
            ret = CPXcallbackgetrelaxationpoint(cb_data, x_p, 0, 0, obj_p)
            obj_p = obj_p[]
            println("Current relaxation point objective: $obj_p.")

            # if (obj_val_non_convex_candidate - obj_p) / obj_val_non_convex_candidate <= max_optimal_physics_gap
            if obj_p > obj_val_non_convex_candidate*1.001
                @info "Prune node because of optimality gap: $obj_p <= $obj_val_non_convex_candidate."
                ret = CPXcallbackprunenode(cb_data)
                return nothing
            end

        end

        # @info "retrieve current solution"
        CPLEX.load_callback_variable_primal(cb_data, context_id)

        π_val =  callback_value.(cb_data, model[:π_avg])
        γ⁺_val = callback_value.(cb_data, model[:γ⁺])
        m⁺_val = callback_value.(cb_data, model[:m⁺])
        γ⁻_val = callback_value.(cb_data, model[:γ⁻])
        m⁻_val = callback_value.(cb_data, model[:m⁻])

        z_val = callback_value.(cb_data, model[:z])
        # Calculate gap and find maximum
        gap = Dict(
            :plus => calc_rel_relaxation_gap.(m⁺_val, π_val, γ⁺_val, Γ_max), # π_val.*γ⁺_val .- m⁺_val.*m⁺_val,
            :minus => calc_rel_relaxation_gap.(m⁻_val, π_val, γ⁻_val, Γ_max)) #π_val.*γ⁻_val .- m⁻_val.*m⁻_val)
        gap_max_side = maximum(gap[:plus]) >= maximum(gap[:minus]) ? Symbol("plus") : Symbol("minus")
        gap_max_index = argmax(Array(gap[gap_max_side]))

        # @info "Gap calculated"
        @info "$(gap[gap_max_side][gap_max_index]) >= $max_relaxation_gap"
        @info "Minimum relaxation gap: $min_max_relaxation_gap."

        # if gap[gap_max_side][gap_max_index] <= min_max_relaxation_gap
        #     min_max_relaxation_gap = gap[gap_max_side][gap_max_index]
        #     @info "Update minimum relaxation gap."
        # end

        println("Best non-convex point objective: $obj_val_non_convex_candidate.")

        if gap[gap_max_side][gap_max_index] >= max_relaxation_gap
            if context_id == CPX_CALLBACKCONTEXT_CANDIDATE
                @info "Reject candidate solution due to relaxation gap."
                ret = CPXcallbackrejectcandidate(
                    cb_data, Cint(0), Cint(0), Cdouble[], Ref{Cchar}(), Ref{Cint}(), Cint[], Cdouble[])
            end

            if context_id == CPX_CALLBACKCONTEXT_BRANCHING

                @info "Side: $gap_max_side, Index: $gap_max_index"
                x_p, obj_p = Vector{Cdouble}(undef, 1), Ref{Cdouble}() # Define empty containeers 
                ret = CPXcallbackgetrelaxationpoint(cb_data, x_p, 0, 0, obj_p)
                obj_p = obj_p[]
                println("Current relaxation point objective: $(obj_p).")

                @info "Start branching."

                column_π = CPLEX.column(cb_data, index(model[:π_avg][gap_max_index]))-1
                column_γ, column_q = gap_max_side == :plus ? 
                    (CPLEX.column(cb_data, index(model[:γ⁺][gap_max_index]))-1, CPLEX.column(cb_data, index(model[:m⁺][gap_max_index]))-1) : 
                    (CPLEX.column(cb_data, index(model[:γ⁻][gap_max_index]))-1, CPLEX.column(cb_data, index(model[:m⁻][gap_max_index]))-1)

                # Select values for specific pipeline and timestep
                γ_val = gap_max_side == :plus ? γ⁺_val[gap_max_index] : γ⁻_val[gap_max_index]
                m_val = gap_max_side == :plus ? m⁺_val[gap_max_index] : m⁻_val[gap_max_index]
                π_val = π_val[gap_max_index]
                # Project (π,γ)-point on to the surface of the cone

                P_sol = [π_val, γ_val, sqrt(π_val*γ_val)]

                # Overwrite solution from last iteration
                # sol_new = (gap_max_side, gap_max_index, [π_val, γ_val, m_val])
                # if sol_new == sol
                #     return println("Pruning: Equivalent branch.")
                # else
                #     sol = sol_new
                # end
 
                @info "P_sol $([π_val, γ_val, m_val])"
                @info "P_sol projected $P_sol"

                # retrieve indices of auxiliaries
                slope = γ_val/π_val
                column_π_b = CPLEX.column(cb_data, index(model[:π_avg_b][gap_max_index]))-1
                column_γ_b = gap_max_side == :plus ? 
                    (CPLEX.column(cb_data, index(model[:γ⁺_b][gap_max_index]))-1) : 
                    (CPLEX.column(cb_data, index(model[:γ⁻_b][gap_max_index]))-1)
    
                π_b_lb, π_b_ub, γ_b_lb, γ_b_ub = 
                    Vector{Cdouble}(undef, 1), Vector{Cdouble}(undef, 1),
                    Vector{Cdouble}(undef, 1), Vector{Cdouble}(undef, 1)
                ret = CPXcallbackgetlocallb(cb_data, π_b_lb, column_π_b, column_π_b)
                ret = CPXcallbackgetlocalub(cb_data, π_b_ub, column_π_b, column_π_b) 
                ret = CPXcallbackgetlocallb(cb_data, γ_b_lb, column_γ_b, column_γ_b) 
                ret = CPXcallbackgetlocalub(cb_data, γ_b_ub, column_γ_b, column_γ_b) 
                π_b_lb, π_b_ub, γ_b_lb, γ_b_ub = π_b_lb[1], π_b_ub[1], γ_b_lb[1], γ_b_ub[1]
    
                @info π_b_lb, π_b_ub, γ_b_lb, γ_b_ub
    
                P_side = [π_b_ub, γ_b_lb, sqrt(π_b_ub*γ_b_lb)]
                coeff_π, coeff_γ, coeff_m = cross(P_side,P_sol) #coeff
    
                @info "Cut L: $([coeff_π, coeff_γ, coeff_m])"
                @info slope
                #### L-Branch
                if slope >= slope_max_point[gap_max_index] #-> shift pi
                    b_var = column_π_b
                    b_val = γ_b_ub/slope # update π_val_lb
                    b_lu = 'L'
                else # -> shift γ
                    b_var = column_γ_b
                    b_val = slope*π_b_ub # update γ_val_ub
                    b_lu = 'U'
                end

                @info b_var, b_val, b_lu

                # Calculate hyperplane
                child_L = Ref{Cint}(0)

                status_L = CPXcallbackmakebranch(
                    cb_data, 1, Cint[b_var], Cchar[b_lu], Cdouble[b_val], 1, 3, Cdouble[cut_slack], Cchar['G'], Cint[0],
                    Cint[column_π, column_γ, column_q], Cdouble[coeff_π, coeff_γ, coeff_m], obj_p[], child_L)

                @info "Added LEFT branch on variable $gap_max_index with node id $(child_L[])."


                P_side = [π_b_lb, γ_b_ub, sqrt(π_b_lb*γ_b_ub)]
                coeff_π, coeff_γ, coeff_m = cross(P_sol,P_side) #coeff
    
                @info "Cut R: $([coeff_π, coeff_γ, coeff_m])"
                #### R-Branch
                if slope >= slope_max_point[gap_max_index] #-> shift pi
                    b_side = column_π_b
                    b_val = γ_b_ub/slope # update π_val_lb
                    b_lu = 'U'
                else # -> shift γ
                    b_side = column_γ_b
                    b_val = slope*π_b_ub # update γ_val_ub
                    b_lu = 'L'
                end

                @info b_var, b_val, b_lu
            
                child_R = Ref{Cint}(0)

                status_R = CPXcallbackmakebranch(
                    cb_data, 1, Cint[b_var], Cchar[b_lu], Cdouble[b_val], 1, 3, Cdouble[cut_slack], Cchar['G'], Cint[0],
                    Cint[column_π, column_γ, column_q], Cdouble[coeff_π, coeff_γ, coeff_m], obj_p[], child_R)


                @info "Add RIGHT branch on variable $gap_max_index with node id $(child_R[])."


                @info "Stop branching."
                if status_L != 0
                    @warn "CPXcallbackmakebranch failed with status $(status_L)"
                elseif status_R != 0
                    @warn "CPXcallbackmakebranch failed with status $(status_U)"
                else
                    @info "Created two new branches.\n The status_U was: $(status_R) (successful). \n The status_L was: $(status_L) (successful)"
                end
            end

        else

            if context_id == CPX_CALLBACKCONTEXT_CANDIDATE # CPX_CALLBACKCONTEXT_BRANCHING


                if gap[gap_max_side][gap_max_index] <= min_max_relaxation_gap
                    min_max_relaxation_gap = gap[gap_max_side][gap_max_index]
                    @info "Update minimum relaxation gap."
                end


                x_p, obj_p = Vector{Cdouble}(undef, 1), Ref{Cdouble}() # Define empty containeers 
                ret = CPXcallbackgetcandidatepoint(cb_data, x_p, 0, 0, obj_p)
                obj_p = obj_p[]

                if obj_p < obj_val_non_convex_candidate
                    @info "Update best bound to $obj_p."
                    obj_val_non_convex_candidate = obj_p

                    # write_solution_sbnb!()

                    γ_val = γ⁺_val .- γ⁻_val
                    m_val = m⁺_val .- m⁻_val
                    @info "z", z_val
                    @info "gap", calc_rel_relaxation_gap.(m_val, π_val, γ_val, Γ_max)
                    @info "gamma", γ_val

                end
            end
    
        end
        return nothing
    end

    MOI.set(model, CPLEX.CallbackFunction(), my_callback_function)
    
    # L_init = Dict(:plus => [[1.,0.,0.] for p in ES.P, t in ES.T],
    #     :minus =>  [[1.,0.,0.] for p in ES.P, t in ES.T])
    # R_init = Dict(:plus => [[0.,1.,0.] for p in ES.P, t in ES.T],
    #     :minus =>  [[0.,1.,0.] for p in ES.P, t in ES.T])
    
    # global sbnb_node_dict = Dict()# Dict{Int64,SbnbNode}() # Collects node information
    # global parent_dict = Dict{Int64,Int64}() # Maps child to parent
    # global last_node = 0

    solve_model!(model)
    write_solution!(ES, model, "MISOCP")

    return nothing
    
end

# ES.q_s = [value(model[:q_s][s,t]) for s in ES.S, t in ES.T]
# ES.q_d_cur = [value(model[:q_d_cur][d,t]) for t in ES.T, d in ES.D]
# ES.m_c = isempty(ES.compressors)==true ? [] : [
#     value(model[:m_c][c,t]) for c in ES.C, t in ES.T]

# ES.m = [value(model[:m][p,t]) for p in ES.P, t in ES.T]




# function compare_bnb_nodes(node1::SbnbNode, node2::SbnbNode)
#     return (node1.p_sol == node2.p_sol) &&
#         (node1.branch_var == node2.branch_var) &&
#         (node1.branch_var_index == node2.branch_var_index)
# end


# function compare_bnb_nodes(node1::Union{SbnbNode,SbnbRootNode}, node2::SbnbRootNode)
#     """
#     Special case for nodes branching from the root node.
#     """
#     return false
# end


# thread_up = 2
# thread_down = 4
# local_progress = 8
# global_progress = 16
# candidate = 32
# relaxation = 64
# branching = 128



# changes = [(CartesianIndex(1,1), [10.,10.,10]), ((CartesianIndex(2,1), [5.,5.,5]))]

# data = deepcopy(L_init[:minus])
# for i in changes
#     data[i[1]] = i[2]
# end

# L_test = update_point(L_init[:minus], CartesianIndex(1,1), [10.,10.,10.])
# sbnb_node_dict[0] = SbnbRootNode(0,L_init,R_init,Dict{Symbol, Float64}())
# parent_dict[0] = 0

# Point (π,γ,m)
