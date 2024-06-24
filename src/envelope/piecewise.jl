function try_piecewise_switches(env::LPEnvelope, r, old_tech, w_limit, A, C_inv; verbose=false)
    
    C_inv .= inv(I(env.n_goods) - (1 + r) * env.A[:, old_tech])

    old_l = copy(vec(env.l[old_tech]))
    w_old = compute_w(C_inv, env.d, old_l)
    w_max = w_old
    l = StrideArray{eltype(C_inv)}(undef, env.n_goods)
    l .= copy(old_l)
    best_sector = best_col = 0

    # These are variables used in compute_w that are preallocated here
    temp_u = StrideArray{eltype(C_inv)}(undef, env.n_goods)
    temp_l_C_inv = StrideArray{eltype(C_inv)}(undef, env.n_goods)
    tech = copy(old_tech)
    
    for sector_tech in eachindex(old_tech) # loop over all sectors
        for country_tech in 1:Int(env.n_countries) # loop over all possible technologies
            # Calculate the column position of the new technology
            new_col = (country_tech - 1) * env.n_goods + sector_tech
            process_old = A[:, old_tech[sector_tech]]
            process_new = A[:, new_col]
            l[sector_tech] = env.l[new_col]
            tech[sector_tech] = new_col

            # Compute the wage resulting from switch to sector_tech
            w = compute_w(C_inv, env.d, l, r, process_old, process_new, sector_tech, temp_u, temp_l_C_inv)
            @inbounds A_trunc = view(env.A, :, tech)
            # w = compute_w(inv(I(36) - (1 + r) * env.A[:, tech]), env.d, l)
            # w = compute_w(env.A[:, tech], I(size(A_trunc, 1)), env.d, Vector(l), r)

            # TODO: Check without Woodbury to account for possible numerical instability
            # TODO: Cache eigval calculation
            # Check whether the wage increased
            if w > w_max
                degenerated = w > w_limit || compute_R(real_eigvals(A_trunc)) < r
                if !degenerated
                    w_max = w
                    best_sector = sector_tech
                    best_col = new_col
                    @debug "Nondegenerated case at profit rate $r"
                else
                    @debug "Degenerated case at profit rate $r"
                end
            end
            # Reset the vectors
            l[sector_tech] = old_l[sector_tech]
            tech[sector_tech] = old_tech[sector_tech]
        end
    end
    if verbose
        println("w_old: $w_old")
        println("w_max: $w_max")
    end
    return w_max > w_old, best_sector, best_col
end