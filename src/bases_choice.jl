
function find_optimal_bases(D::Number=0.5, ΔT::Number=3600, ΔX::Number=5000, c::Number=0)
    m = Model(Ipopt.Optimizer)
    @variables m begin
        D_base >= 0
        T_base >= 0
    end
    @NLobjective(m, Min,
        (log10(T_base/ΔT)-c)^2 + 
        (log10((350^2*T_base^2*D_base)/ΔX)-c)^2 +
        (log10(T_base/ΔT)-c)^2 +
        (log10(1/(D_base*ΔX))-c)^2 +
        (log10(0.01*350^2*T_base^2*D_base/(2*D^3*1))-c)^2
    )
    optimize!(m)

    return value(D_base), value(T_base)
end

# D, T = find_optimal_bases(0.5, 3600, 5000, 0)