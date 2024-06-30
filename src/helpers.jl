""" Compute the real valued eigenvalues. """
real_eigvals(x) = real.(filter(e -> imag(e) == 0, eigvals(x)))

"""
Compute the highest sensible profit rate for given eigenvalues of the input-output matrix.

Briefly stated, an application of the Perron-Frobenius theorem gives the existence of a
maximum profit rate R_max which gives non-negative prices between [0, R_max].
"""
function compute_R(eigenvalues)
    eig_frob = maximum(eigenvalues)
    return (1.0 - eig_frob) / eig_frob
end

""" For normalization it is conventient to return 0.0 for division by zero. """
safe_divide(x, y) = y == 0 ? 0.0 : x / y

replace_with_zero(x, value = 1.0) = x == value ? 0.0 : x

""" Setup and return a solver for our linear problem with intensities that can modify the profit rate. """
function create_intensities_r_solver(solver, l, A, B, d, lb)
    model = direct_model(solver)
    @variable(model, x[i = eachindex(lb)]>=lb[i])
    @variable(model, y[i = eachindex(lb)])
    @objective(model, Min, l*x)
    @constraint(model, con[i = eachindex(lb)], 1 * x[i] - y[i]==0)
    @constraint(model, d, B * x - A * y.>=d)
    return model
end

""" Setup and return a solver for our linear problem with intensities. """
function create_intensities_solver(solver, l, C, d, lb)
    model = direct_model(solver)
    @variable(model, x[i = eachindex(lb)].>=lb[i])
    @objective(model, Min, sum(x[i] * l[i] for i in eachindex(lb)))
    @constraint(model, con, C * x.≥d)
    return model
end

""" Setup and return a solver for our linear problem with prices. """
function create_prices_solver(solver, l, C, d, lb)
    model = direct_model(solver)
    @variable(model, p[i = eachindex(lb)].>=lb[i])
    @objective(model, Max, sum(p[i] * d[i] for i in eachindex(lb)))
    @constraint(model, con, C' * p.≤vec(l))
    return model
end

""" Modify the input-matrix to A. """
function modify_A!(model, A, B, d)
    delete(model, model[:d])
    unregister(model, :d)
    @constraint(model, d, B * model[:x] - A * model[:y].>=d)
end

""" Set the coefficient of variable[j] in the i-th constraint to C[i, j]. """
function modify_C!(model, variable, C)
    for j in axes(C, 2) # Loop over the columns of C (each variable)
        set_normalized_coefficient.(model[:con], variable[j], vec(C[:, j]))
    end
end

""" Set a new profit rate r. """
function modify_r!(model, r)
    set_normalized_coefficient(model[:con], model[:x], fill(1 + r, length(model[:x])))
end

""" Calculate the wage rate w.

Since A and l are transposed compared to the usual formula (A from the start, l inside this function), we need to reverse the order.
"""
compute_w(A, B, d, l, r) = 1 / (l' * ((B - (1 + r) * A) \ d))

""" Calculate the wage rate w for a given inverse matrix.
"""
compute_w(C_inv, d, l) = 1 / (l' * C_inv * d)

""" Calculate the wage rate w for a given LU decomposition.
"""
compute_w(d, l; lu) = 1 / (l' * (lu \ d))

""" Calculate the wage rate w with the Woodbury formula.

l is a column vector that gets transposed.
"""
function compute_w(C_inv, d, l, r, process_old, process, industry, u, l_C_inv)
    # We substract -(1+r)A from B, therefore we have to reverse the sign here
    @.. u = (process_old - process) * (1 + r)
    # This is v' * C_inv given v = eachindex(d) .== industry
    v_T_C_inv = view(C_inv, industry, :)
    # It would be smarter to look at C_inv * d since this remains constant
    mul!(l_C_inv, C_inv', l)
    denom = 1 + dot(v_T_C_inv, u)
    numerator = dot(l_C_inv, u) * dot(v_T_C_inv, d)
    return 1 / (dot(l_C_inv, d) - numerator / denom) # dot is needed bc of StrideArrays
end

""" Calculate the wage rate w given a LU decomposition using the Woodbury formula.

l is a column vector that gets transposed.

Plugin Woodbury: (A + uv')^{-1} d =: x_new
u_new = A^{-1} u
This gives: x_new = x_old - u_new v^T x_old / (1 + v^T u_new)
Result: 1 / (x_new^T * d)

See also https://scicomp.stackexchange.com/questions/21303/solve-rank-one-update-to-lu-using-plain-vanilla-lu-routine
"""
function compute_w(decomp, d, l, r, process_old, process, industry, u, d_C_inv_old, u_new)
    # We substract -(1+r)A from B, therefore we have to reverse the sign here
    @.. u = (process_old - process) * (1 + r)
    ldiv!(u_new, decomp, u)
    denom = 1 + u_new[industry]
    # Straightforward implementation:
    # numerator = u_new * d_C_inv_old[industry]
    # return 1 / (dot(l, d_C_inv_old - numerator / denom)) # dot is needed bc of StrideArrays
    s = 0.0
    for i in eachindex(l)
        woodbury = d_C_inv_old[i] - u_new[i] * d_C_inv_old[industry] / denom
        s += woodbury * l[i]
    end
    return 1 / s
end

# TODO: Think about utilizing the intensities to get some kind if derivative wrt. direction of new process.