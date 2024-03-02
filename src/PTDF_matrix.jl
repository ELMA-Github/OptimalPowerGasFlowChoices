function create_PTDF_matrix(L_df)
    # determination of dimensions requires "proper" data matrix
    no_nodes = maximum(vcat(L_df[!, :Start], L_df[!, :Stop]))
    no_lines = maximum(L_df[!,:Line_num])
    # Incidence matrix C
    C = zeros(no_nodes, no_lines)
    for row in 1:nrow(L_df)
        C[L_df[row, :Start], row] = 1
        C[L_df[row, :Stop], row] = -1
    end
    # Bus susceptance matrix B
    B = diagm(1 ./ L_df[!,:X_pu]) # susceptance matrix
    L = C*B*C' # Laplacian matrix

    Lbar = L[1:end-1, 1:end-1]
    Btilde = zeros(no_nodes,no_nodes)
    Btilde[1:end-1, 1:end-1] = inv(Lbar)
    
    return B*C'*Btilde

end