global precision_r_diff = 1e-6

# TODO: ensure that the techs have all the same length
# TODO: use B properly if not identity
struct BinaryEnvelope{T, R}
    A
    B
    l
    d
    n_goods::Int64
    n_countries::Int64
    ids::Dict{Int64, Vector{T}} # map technologies to ids
    tech_dict::Dict{R, Int64} # map profit rates to ids
    reverse::Dict{Int64, Set{R}} # map ids to profit rates
    wages::Dict{R, Float64} # map to wages
    stepsize::Float64
    Envelope(; A, B, l, d, profit_rates, initial_tech, stepsize=0.0001) = new{Int64, eltype(profit_rates)}(
        A,
        B,
        l,
        d,
        size(A, 1),
        div(size(A, 2), n_goods),
        Dict(1 => initial_tech),
        Dict(profit_rates .=> 1),
        Dict(1 => Set(profit_rates)),
        Dict(profit_rates .=> -Inf),
        stepsize
    )
end

# Update the envelope with a new technology
function update_envelope(envelope, r, technology, w)
    # set the new technology
    old_id = envelope.tech_dict[r]
    new_id = findfirst(==(technology), envelope.ids)
    @assert old_id != new_id "You have to update with a different technology."
    println(new_id)
    if isnothing(new_id)
        new_id = length(envelope.ids) + 1
        envelope.ids[new_id] = technology
        envelope.reverse[new_id] = Set([r])
    end
    envelope.tech_dict[r] = new_id
    push!(envelope.reverse[new_id], r)
    # Remove the profit rate from the reverse dict
    println(old_id, new_id)
    pop!(envelope.reverse[old_id], r)
    # set the new maximal wage
    envelope.wages[r] = w
end

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

function try_piecewise_switches(envelope::LPEnvelope, r, C_inv, old_tech)
    old_l = vec(envelope.l[old_tech])
    w_old = compute_w(C_inv, envelope.d, old_l)
    w_max = w_old
    l = copy(old_l)
    best_sector = best_col = 0

    # These are variables used in compute_w that are preallocated here
    temp_u = Vector{Float64}(undef, envelope.n_goods)
    temp_l_C_inv = adjoint(Vector{Float64}(undef, envelope.n_goods))
    
    for sector_tech in eachindex(old_tech) # loop over all sectors
        for country_tech in 1:Int(envelope.n_countries) # loop over all possible technologies
            # Calculate the column position of the new technology
            new_col = (country_tech - 1) * envelope.n_goods + sector_tech
            process_old = view(envelope.A, :, old_tech[sector_tech])
            process_new = view(envelope.A, :, new_col)
            l[sector_tech] = envelope.l[new_col]
            # Compute the wage resulting from switch to sector_tech
            w = compute_w(C_inv, envelope.d, l, r, process_old, process_new, sector_tech, temp_u, temp_l_C_inv)

            # Check whether the wage increased
            if w > w_max
                w_max = w
                best_sector = sector_tech
                best_col = new_col
            end
            # Reset the l vector
            l[sector_tech] = old_l[sector_tech]
        end
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

function binary_search(envelope, start_r, end_r, start_terminated=false)
    if end_r - start_r < precision_r_diff
        # TODO:
        # check for numerical instability
    else
        # First, we try to improve the current tech at start_r (if not already optimal)
        # If tech transfer works:
        #   If piecewise improvement in end_r is possible: Reverse roles.
        #   If not: Check all values between start_r and end_r (if necessary, search in smaller interval again)
        # If it doesn't work: split in lower and upper half of the profit rates and try again
        tech = envelope.ids[envelope.tech_dict[start_r]]
        if !start_terminated
            C_inv = inv(B[:, tech] - (1+r) * A[:, tech])
            start_terminated = try_piecewise_switches(envelope, start_r, d, C_inv)
        end
        # Try to use newly found tech and update the inverse
        transfer_worked, C_inv_new = try_start_tech(envelope, start_r, end_r)
        if transfer_worked
            end_terminated = try_piecewise_switches(envelope, end_r, d, C_inv_new)
            all_terminated = end_terminated
            # If we are done at this profit rate: Check transfer to next profit rates and if they are done as well
            if end_terminated
                next_r = r_end + envelope.stepsize
                for interval_r in next_r:flipsign(envelope.stepsize, next_r - r_start):r_start
                    transfer_worked, C_inv_interv = try_start_tech(envelope, end_r, interval_r)
                    also_terminated = try_piecewise_switches(envelope, interval_r, d, C_inv_interv)
                    # If we are not done: search in this smaller interval again
                    if !also_terminated
                        start_r = interval_r
                        end_r = start_r
                        all_terminated = false
                        break
                    end
                end                
            end
            all_terminated && return
            # Maybe here extend if end_r is highest, but this is tricky
            binary_search(envelope, C_inv_new, end_r, start_r)
        else
            # Maybe rationals/fixed point arithmetic would be better
            # Split the interval in half
            grid = start_r:envelope.precision:end_r
            r = grid[div(length(grid), 2)]

            # Try the lower half (without new pairwise search in start_r if we are already optimal)
            binary_search(envelope, start_r, r, start_terminated)
            # Try the upper half
            binary_search(envelope, end_r, r + envelope.stepsize)
        end
    end
end

function test_at_corners(envelope, l, d, R)
    # test at the lowest profit rate
    # try_piecewise_switches(envelope, 0.0, l, d, ...)
    # test at the highest profit rate
    # try_piecewise_switches(envelope, R, l, d, ...)
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