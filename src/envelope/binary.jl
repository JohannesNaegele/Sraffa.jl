function try_piecewise_switches(envelope::BinaryEnvelope, r, C_inv)
    w_old = envelope.wages[r]
    w_max = w_old
    old_tech_id = envelope.tech_dict[r]
    old_l = envelope.l[envelope.ids[old_tech_id]]
    l = copy(old_l)
    for sector_tech in 1:envelope.n_goods
        for country_tech in 1:envelope.n_countries
            process_old = view(envelope.A, :, envelope.ids[old_tech_id][sector_tech])
            new_col = (country_tech - 1) * envelope.n_goods + sector_tech
            process_new = view(envelope.A, :, new_col)
            l[sector_tech] = envelope.l[new_col]
            w = compute_w(C_inv=C_inv, d=envelope.d, l=l, process_old=process_old, process=process_new, industry=sector_tech, r=r)
            w > w_max && (w_max = w)
        end
        l[sector_tech] = old_l[sector_tech] # reset
    end
    # Update the envelope
    update_envelope(envelope, r, new_tech, w_max)
    return w_max != w_old
end

function try_piecewise_switches(env::LPEnvelope, r, old_tech, w_limit, A, C_inv; verbose=false)
    
    C_inv .= inv(I(36) - (1 + r) * env.A[:, old_tech])

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

            # Compute the wage resulting from switch to sector_tech
            w = compute_w(C_inv, env.d, l, r, process_old, process_new, sector_tech, temp_u, temp_l_C_inv)
            tech[sector_tech] = new_col
            @inbounds A_trunc = view(env.A, :, tech)
            # w = compute_w(inv(I(36) - (1 + r) * env.A[:, tech]), env.d, l)
            # w = compute_w(A, I(36), env.d, l, r)

            # Check whether the wage increased
            if w > w_max && w < w_limit # TODO: does this work?
                if !(compute_R(real_eigvals(A_trunc)) < r)
                    w_max = w
                    best_sector = sector_tech
                    best_col = new_col
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

function try_start_tech(envelope, start_r, end_r)
    tech = envelope.ids[envelope.tech_dict[start_r]]
    C_inv = envelope.B - (1 + end_r) * envelope.A[:, tech]
    w_old = envelope.wages[end_r]
    process_old = view(envelope.A, :, envelope.ids[old_tech_id][sector_tech])
    new_col = (country_tech - 1) * n_goods + sector_tech
    process_new = view(envelope.A, :, new_col)
    l[sector_tech] = envelope.l[new_col]
    # It should be ok to update the inverse with the Woodbury formula since the new process is usually quite different
    w, C_inv_new = compute_w(C_inv=C_inv, d=envelope.d, l=l, process_old=process_old, process=process_new, industry=sector_tech, r=end_r - start_r)
    update_envelope(envelope, r, new_tech, w_max)
    return w_max != w_old, C_inv_new
end

""" Compute the envelope for the book-of-blueprints A with the VFZ algorithm.

It might make sense to at least sometimes use the Sherman-Morrison formula to update the inverse.

This left-right stuff just means that we don't want to try out again all possibilities at the next higher/lower step but use the new best from above/below.
Left/right side are sufficient for induction since when they land on the envelope then we have piecewise switches guaranteed.
That is also the reason why numerical instability will not be a problem since the corners are checked piecewise.
The idea is that after left-to-right on the right we have better tech, therefore we go in the opposite direction.

My idea is to use something along the lines of binary search. We start with the lowest profit rate and then we try to find the best wage for piecewise switches.
Then we try this technology at the highest profit rate which maps to the same technology. If it doesn't work we try all possibilities there and then proceed to the profit rate in the middle.
For that we can again use Sherman-Morrison to update the inverse. This should also be numerically stable for strongly varying profit rates.
At the end we should check again at lowest and highest profit rate of one wage curve on the envelope for numerical instability and redo binary without Woodbury.

Bear in mind that the columns of A describe the technology of one sector.
Therefore all objects have permuted dimensions compared to the usual ``(1 + r)Ap + wl = p`` formula!!
"""
function compute_vfz(; A, B, l, d, R, step, verbose = false, save_all=true)
    # Number of goods
    n_goods = size(A, 1)
    profit_rates = 0.0:step:R
    # Initialize with technology of first country
    init_tech = A[:, 1:n_goods]

    envelope = Envelope(
        A=A,
        B=B,
        l=l,
        d=d,
        profit_rates=profit_rates,
        initial_tech=init_tech
    )

    binary_search(envelope, 0.0, R)
    # Extend the maximal profit rate
    # compute_R(maximum.(real_eigvals.([A1, A2])))
    # Test whether the Woodbury formula is stable; do this only at begin/start of tech choice
    test_at_corners(envelope, l, d, R)
    n_techs = unique(values(envelope.tech_dict))
    println("We have $(n_techs - 1) switches.")
    return envelope
end