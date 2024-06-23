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