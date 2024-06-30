""" Compute the envelope for the book-of-blueprints A with linear programming.

Bear in mind that the columns of A describe the technology of one sector.
Therefore all objects have permuted dimensions compared to the usual ``(1 + r)Ap + wl = p`` formula!!

We use all sectors for switchpoint calculation but consider only `effects_sectors` for effects such as reverse capital deepening.
"""
function compute_envelope(; A, B, l, d, R, stepsize, model_intensities,
        model_intensities_trunc, model_prices, save_all=true, effects_sectors=1:size(A, 1), extend=true, calc_all=true)

    @assert all(d .>= 0) "d must be non-negative!"

    r_max_old = -Inf
    r_max_new = Inf
    i_old = 1

    # Grid of profit rates
    profit_rates = 0:stepsize:R

    env = LPEnvelope(A, l, d, profit_rates, effects_sectors)

    C_inv_temp = similar(A[:, 1:size(A, 1)])
    A = StrideArray{eltype(env.A)}(undef, size(env.A))
    copyto!(A, env.A)

    while r_max_old < r_max_new
        for (i, r) in enumerate(profit_rates)
            # Get the indices for the columns with the highest intensity per sector
            function highest_intensity_indices(sector)
                max_index = argmax(j -> env.intensities[j, i], sector:env.n_goods:length(l))
                return max_index
            end
            if r > r_max_old
                @debug "Profit rate: $r"
                # TODO: Regularly check with LP case whether piecewise switches succeeded
                if calc_all || i == 1
                    # Compute intensities
                    modify_r!(model_intensities, r)
                    optimize!(model_intensities) # min l * x s.t. C * x ≥ d AND x .≥ 0
                    error_msg = "$(termination_status(model_intensities))\nProbably unfeasible d!"
                    @assert is_solved_and_feasible(model_intensities) error_msg
                    env.intensities[:, i] = replace_with_zero.(value.(model_intensities[:x]))
                    map!(sector -> highest_intensity_indices(sector), view(env.chosen_technology, :, i), 1:env.n_goods)
                    env.chosen_technology[:, i] = Main.zamb_res[:, 1]
                else
                    # For pairwise switches use last optimal tech
                    tech = view(env.chosen_technology, :, i - 1)
                    copyto!(view(env.chosen_technology, :, i), tech)
                    w_limit = compute_w(A[:, tech], B[:, tech], d, l[tech], profit_rates[i - 1])
                    improved, best_sector, best_col = try_piecewise_switches(env, r, tech, w_limit, A, C_inv_temp)
                    if improved
                        env.chosen_technology[best_sector, i] = best_col
                    end
                end
            end
        end

        r_max_old = profit_rates[end]
        # FIXME: This only works for B being the identity matrix
        r_max_new = compute_R(real_eigvals(A[:, env.chosen_technology[:, end]]))
        i_old = length(profit_rates)
        if extend
            profit_rates = 0:stepsize:r_max_new
        end
        # This is effectively rounding to ensure that old equals new if we are advancing sub-stepsize
        r_max_new = profit_rates[end]
        if r_max_old < r_max_new
            extend!(env, length(profit_rates) - i_old)
        end
    end

    d = ones(env.n_goods)
    set_normalized_rhs.(model_intensities_trunc[:con], d)
    set_objective_coefficient.(model_prices, model_prices[:p], d)

    # We compute the truncated stuff and prices only at switch points
    for (i, tech) in enumerate(eachcol(env.chosen_technology)[begin:(end - 1)])
        # Check whether we are at switch point
        switch = tech[effects_sectors] != env.chosen_technology[effects_sectors, i + 1]
        if switch
            for j in i:(i + 1)
                # FIXME: Do we want this at all?
                # Compute truncated intensities
                A_trunc = view(env.A, :, env.chosen_technology[:, j])  # filter columns
                l_trunc = vec(view(l, :, env.chosen_technology[:, j]))'  # filter columns
                B_trunc = I(size(A_trunc, 1))  # initialize identity matrix
                C_trunc = B_trunc - A_trunc
                modify_C!(model_intensities_trunc, model_intensities_trunc[:x], C_trunc)
                optimize!(model_intensities_trunc)  # min l_trunc * x s.t. C_trunc * x ≥ d AND x .≥ 0
                @assert is_solved_and_feasible(model_intensities_trunc) "switch at r=$(profit_rates[i])"
                env.intensities_trunc[:, j] = value.(model_intensities_trunc[:x])

                # Compute prices
                r = profit_rates[j]
                C = (B - (1 + r) * A)
                C_trunc = view(C, :, env.chosen_technology[:, j])  # filter columns
                set_normalized_rhs.(model_prices[:con], vec(l_trunc))
                modify_C!(model_prices, model_prices[:p], permutedims(C_trunc))
                optimize!(model_prices)  # max p * d s.t. C' * p ≤ l AND p .≥ 0
                @assert is_solved_and_feasible(model_prices)
                env.prices[:, j] = value.(model_prices[:p])

                # Compute (l * x)
                env.lx[j] = l_trunc[:, effects_sectors] * env.intensities_trunc[:, j][effects_sectors]
                # Compute (pA)'
                env.pA[:, j] = vec(A_trunc[effects_sectors, effects_sectors]' * env.prices[:, i][effects_sectors])
                # Compute A * x / (l * x)
                env.right_side_factor[:, j] = A_trunc[effects_sectors, effects_sectors] *
                    env.intensities_trunc[:, j][effects_sectors] / env.lx[j]
            end

            # Compute p * A * x / (l * x)
            push!(env.capital_intensities,
                [
                    profit_rates[i] => env.prices[:, i][effects_sectors]' * env.right_side_factor[:, i],
                    profit_rates[i + 1] => env.prices[:, i][effects_sectors]' * env.right_side_factor[:, i + 1]
                ]
            )
            push!(env.lx_at_switch, [
                    profit_rates[i] => env.lx[i],
                    profit_rates[i + 1] => env.lx[i + 1]
                ]
            )
            push!(env.pA_at_switch,
                [
                    profit_rates[i] => env.pA[:, i][effects_sectors],
                    profit_rates[i + 1] => env.pA[:, i + 1][effects_sectors]
                ]
            )
            push!(env.l_at_switch,
                [
                    profit_rates[i] => vec(view(l, :, env.chosen_technology[:, i]))[effects_sectors],
                    profit_rates[i + 1] => vec(view(l, :, env.chosen_technology[:, i + 1]))[effects_sectors]
                ]
            )
            push!(env.intensities_at_switch,
            [
                profit_rates[i] => env.intensities_trunc[:, i][effects_sectors],
                profit_rates[i + 1] => env.intensities_trunc[:, i + 1][effects_sectors]
            ]
            )
            push!(env.prices_switch, profit_rates[i] => env.prices[:, i][effects_sectors])
            push!(env.technologies_switch,
                [profit_rates[i] => env.chosen_technology[:, i][effects_sectors],
                profit_rates[i + 1] => env.chosen_technology[:, i + 1][effects_sectors]]
            )
        end
    end

    profit_rates = 0:stepsize:r_max_new

    # Precision used to round to string
    precision = Int(ceil(log10(1 / stepsize)))

    profit_rates_to_names = Dict(profit_rates .=>
        string.(round.(profit_rates, digits = precision)))
    df_intensities = DataFrame(
        env.intensities, map(r -> profit_rates_to_names[r], profit_rates))

    switches = save(env, save_all)

    return df_intensities, profit_rates_to_names, profit_rates, switches
end