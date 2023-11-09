function save_results_to_excel(ES, config_dict)
    # Save to EXCEL
    if isempty(config_dict[:warmstart])
        warm ="NO"
    else
        warm = config_dict[:warmstart][:method]
    end

    if config_dict[:space_disc] == 0
        dx ="pipe"
    else
        dx =config_dict[:segment]
    end

    name_file=string("output_",config_dict[:case_study],"_",config_dict[:model_type],"_PDE", config_dict_algorithm[:PDE_type],"_OB",config_dict[:objective_type],"_warm",warm,"_dt",config_dict[:dt],"_dx",dx,"_pu",config_dict[:pu],".xlsx")
    output_path=joinpath(dirname(@__DIR__), "results", name_file)

    # Optimization time and objective value
    gas_cost_ex_post = ES.bases_dict[:Cost_base]* total_cost_gas(ES.objective_type, ES.dt, ES.q_s, ES.sources, ES.q_d_cur, ES.C_cur_load_gas, ES.π, ES.compressors, ES.C, ES.S, ES.D_gas, ES.T)
    df = DataFrame(obj_val = [gas_cost_ex_post], comp_time = [ES.solve_time]) #works only if cost_base= 1 or if linear ob value
    # Pressure
    df_pressure = DataFrame(ES.π*ES.bases_dict[:Π_base],["t$i" for i in 1:size(ES.π, 2)])
    df_pressure  = insertcols!(df_pressure , 1, :Node =>["n$i" for i in 1:size(ES.π, 1)])
    #XLSX.writetable(output_path, df_pressure ,  sheetname="Pressure",overwrite=true)
    # Gas source
    df_source = DataFrame(ES.q_s*ES.bases_dict[:M_base],["t$i" for i in 1:size(ES.q_s, 2)])
    df_source  = insertcols!(df_source , 1, :Source =>["s$i" for i in 1:size(ES.q_s, 1)])
    # Gas flow
    df_flow = DataFrame(ES.m*ES.bases_dict[:M_base],["t$i" for i in 1:size(ES.m, 2)])
    df_flow  = insertcols!(df_flow , 1, :Pipe =>["p$i" for i in 1:size(ES.m, 1)])

    # Gas flow in
    df_flow_in = DataFrame(ES.m_in*ES.bases_dict[:M_base],["t$i" for i in 1:size(ES.m_in, 2)])
    df_flow_in  = insertcols!(df_flow_in , 1, :Pipe =>["p$i" for i in 1:size(ES.m_in, 1)])

    # Gas flow out
    df_flow_out = DataFrame(ES.m_out*ES.bases_dict[:M_base],["t$i" for i in 1:size(ES.m_out, 2)])
    df_flow_out  = insertcols!(df_flow_out , 1, :Pipe =>["p$i" for i in 1:size(ES.m_out, 1)])

    # Gas compressor
    if ~isempty(ES.m_c)
        df_compr = DataFrame(ES.m_c*ES.bases_dict[:M_base],["t$i" for i in 1:size(ES.m_c, 2)])
        df_compr  = insertcols!(df_compr , 1, :Compressor =>["c$i" for i in 1:size(ES.m_c, 1)])
    else
        df_compr = DataFrame([Vector{Union{Missing, Float64}}(missing, 0) for _ in 1:ES.N_dt], ["t$i" for i in 1:ES.N_dt])
    end

    # Curtailment
    df_curt = DataFrame(ES.q_d_cur*ES.bases_dict[:M_base],["t$i" for i in 1:size(ES.q_d_cur, 2)])
    df_curt  = insertcols!(df_curt , 1, :Demand =>["d$i" for i in 1:size(ES.q_d_cur, 1)])

    # Gas Load
    df_gasload = DataFrame([Vector{Union{Missing, Float64}}(missing, 0) for _ in 1:ES.N_dt], ["t$i" for i in 1:ES.N_dt])
    for d_gas in ES.D_gas
        push!(df_gasload,ES.demands[d_gas].Q_d[1:ES.N_dt]*ES.bases_dict[:M_base])
    end
    df_gasload = insertcols!(df_gasload, 1, :Demand =>["d$i" for i in 1:size(ES.D_gas, 1)])

    # Linepack
    lp = calc_linepack_expost(ES)*ES.bases_dict[:Π_base]*ES.bases_dict[:D_base]*ES.bases_dict[:T_base]^2*10^6;

    df_lp = DataFrame(lp,["t$i" for i in 1:size(lp, 2)])
    df_lp  = insertcols!(df_lp, 1, :Pipe =>["p$i" for i in 1:size(lp, 1)])

    # Linepack difference

    lp_diff = diff(lp,dims=2);
    df_lp_diff = DataFrame(lp_diff,["t$i" for i in 2:size(lp, 2)])
    df_lp_diff  = insertcols!(df_lp_diff, 1, :Pipe =>["p$i" for i in 1:size(lp_diff, 1)])

    # Relative relaxation gap
    gap_rel = calc_rel_relaxation_gap(ES);
    df_gap_rel  = DataFrame(gap_rel,["t$i" for i in 1:size(gap_rel, 2)])
    df_gap_rel  = insertcols!(df_gap_rel , 1, :Pipe =>["p$i" for i in 1:size(gap_rel, 1)])

    if sys == "gas_only"


        XLSX.writetable(output_path, overwrite=true,
            OB_CompTime=(collect(eachcol(df)), names(df)),
            Pressure=(collect(eachcol(df_pressure)), names(df_pressure)),
            GasSource=(collect(eachcol(df_source)), names(df_source)),
            GasFlow=(collect(eachcol(df_flow)), names(df_flow)),
            GasFlowIn=(collect(eachcol(df_flow_in)), names(df_flow_in)),
            GasFlowOut=(collect(eachcol(df_flow_out)), names(df_flow_out)),
            GasCompressor=(collect(eachcol(df_compr)), names(df_compr)),
            GasCurtailment=(collect(eachcol(df_curt)), names(df_curt)),
            GasLoad=(collect(eachcol(df_gasload)), names(df_gasload)),
            Linepack=(collect(eachcol(df_lp)), names(df_lp)),
            LinepackDiff=(collect(eachcol(df_lp_diff)), names(df_lp_diff)),
            RelaxationRel= (collect(eachcol(df_gap_rel )), names(df_gap_rel )), 
        )

    elseif sys== "integrated_power_and_gas"

        # Optimization time and objective value
        gas_cost_ex_post = ES.bases_dict[:Cost_base]* total_cost_gas(ES.objective_type, ES.dt, ES.q_s, ES.sources, ES.q_d_cur, ES.C_cur_load_gas, ES.π, ES.compressors, ES.C, ES.S, ES.D_gas, ES.T)
        obj_power_ex_post = ES.power_cost*ES.bases_dict[:Cost_base]
        df_int = DataFrame(obj_val = [gas_cost_ex_post+obj_power_ex_post], obj_gas = [gas_cost_ex_post], obj_power = [obj_power_ex_post], comp_time = [ES.solve_time])

        # Gas Generator consumption
        df_gasgen = DataFrame(ES.q_g*ES.bases_dict[:M_base],["t$i" for i in 1:size(ES.q_g, 2)])
        df_gasgen  = insertcols!(df_gasgen , 1, :GasGenerator =>["gfpp$i" for i in 1:size(ES.q_g, 1)])

        # Generator
        df_gen = DataFrame(ES.p*ES.bases_dict[:S_base],["t$i" for i in 1:size(ES.p, 2)])
        df_gen  = insertcols!(df_gen , 1, :GasGenerator =>["g$i" for i in 1:size(ES.p, 1)])

        # Wind dispatch
        df_wind = DataFrame(ES.w*ES.bases_dict[:S_base],["t$i" for i in 1:size(ES.w, 2)])
        df_wind  = insertcols!(df_wind , 1, :Wind =>["w$i" for i in 1:size(ES.w, 1)])

        # Wind production
        df_wind_prod = DataFrame([Vector{Union{Missing, Float64}}(missing, 0) for _ in 1:ES.N_dt], ["t$i" for i in 1:ES.N_dt])
        for wind in ES.J
            push!(df_wind_prod,ES.windgenerators[wind].P_w[1:ES.N_dt]*ES.bases_dict[:S_base])
        end
        df_wind_prod = insertcols!(df_wind_prod, 1, :Wind =>["w$i" for i in 1:size(ES.J, 1)])

        # Power demand
        df_powerload = DataFrame([Vector{Union{Missing, Float64}}(missing, 0) for _ in 1:ES.N_dt], ["t$i" for i in 1:ES.N_dt])
        for d in ES.D_el
            push!(df_powerload,ES.demands_EL[d].D_el[1:ES.N_dt]*ES.bases_dict[:S_base])
        end
        df_powerload  = insertcols!(df_powerload , 1, :Demand =>["d$i" for i in 1:size(ES.D_el, 1)])

        # Power Curtailment
        df_powercur = DataFrame(ES.p_d_cur*ES.bases_dict[:S_base],["t$i" for i in 1:size(ES.p_d_cur, 2)])
        df_powercur = insertcols!(df_powercur , 1, :Wind =>["d$i" for i in 1:size(ES.p_d_cur, 1)])

        XLSX.writetable(output_path, overwrite=true,
            OB_CompTime=(collect(eachcol(df_int)), names(df_int)),
            Pressure=(collect(eachcol(df_pressure)), names(df_pressure)),
            GasSource=(collect(eachcol(df_source)), names(df_source)),
            GasFlow=(collect(eachcol(df_flow)), names(df_flow)),
            GasFlowIn=(collect(eachcol(df_flow_in)), names(df_flow_in)),
            GasFlowOut=(collect(eachcol(df_flow_out)), names(df_flow_out)),
            GasCompressor=(collect(eachcol(df_compr)), names(df_compr)),
            GasCurtailment=(collect(eachcol(df_curt)), names(df_curt)),
            GasLoad=(collect(eachcol(df_gasload)), names(df_gasload)),
            Linepack=(collect(eachcol(df_lp)), names(df_lp)),
            LinepackDiff=(collect(eachcol(df_lp_diff)), names(df_lp_diff)),
            RelaxationRel= (collect(eachcol(df_gap_rel )), names(df_gap_rel )), 
            GasGenerator=(collect(eachcol(df_gasgen)), names(df_gasgen)),
            Generator=(collect(eachcol(df_gen)), names(df_gen)),
            WindDispatch = (collect(eachcol(df_wind)), names(df_wind)),
            WindProduction = (collect(eachcol(df_wind_prod)), names(df_wind_prod)),
            PowerLoad = (collect(eachcol(df_powerload)), names(df_powerload)),
            PowerCurtailment = (collect(eachcol(df_powercur)), names(df_powercur))
    )

    end
    return nothing
end