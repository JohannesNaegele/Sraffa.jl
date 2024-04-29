Base.@kwdef struct LPEnvelope <: EnvelopeMethod
    # Number of goods
    n_goods

    # Activity levels (q)
    intensities
    # Truncated activity levels; calculated via C = (B - A)
    intensities_trunc
    pA
    right_side_factor
    lx
    # Prices computed from truncated intensities via the dual representation
    prices
    # Indices of chosen technologies; chosen_technology[i] ∈ i * n_goods
    chosen_technology

    # These are data structures used to store only at switches
    # Capital intensities (pAx / lx)
    capital_intensities
    lx_at_switch
    intensities_at_switch
    # Prices (specified to be exactly the same before/after switch; it doesn't make a difference though)
    prices_switch
    technologies_switch
    pA_at_switch
    l_at_switch
end


function LPEnvelope(A, profit_rates, effects_sectors)
    n_goods = size(A, 1)
    LPEnvelope(
        n_goods = n_goods,
        intensities = zeros(size(A, 2), length(profit_rates)),
        intensities_trunc = zeros(size(A, 1), length(profit_rates)),
        pA = zeros(effects_sectors, length(profit_rates)),
        right_side_factor = zeros(effects_sectors, length(profit_rates)),
        lx = zeros(length(profit_rates)),
        prices = zeros(n_goods, length(profit_rates)),
        chosen_technology = zeros(Int64, n_goods, length(profit_rates)),
        capital_intensities = Vector{Pair{Float64, Float64}}[],
        lx_at_switch = Vector{Pair{Float64, Float64}}[],
        intensities_at_switch = Vector{Pair{Float64, Vector{Float64}}}[],
        prices_switch = Pair{Float64, Vector{Float64}}[],
        technologies_switch = Vector{Pair{Float64, Vector{Int64}}}[],
        pA_at_switch = Vector{Pair{Float64, Vector{Float64}}}[],
        l_at_switch = Vector{Pair{Float64, Vector{Float64}}}[]
    )
end

""" Compute the envelope for the book-of-blueprints A with linear programming.

Bear in mind that the columns of A describe the technology of one sector.
Therefore all objects have permuted dimensions compared to the usual ``(1 + r)Ap + wl = p`` formula!!

We use all sectors for switchpoint calculation but consider only `effects_sectors` for effects such as reverse capital deepening.
"""
function compute_envelope(; A, B, l, d, R, step, model_intensities,
        model_intensities_trunc, model_prices, verbose = false, save_all=true, effects_sectors=33)

    # Precision used to round to string
    precision = Int(ceil(log10(1 / step)))

    R_max = Inf
    # Grid of profit rates
    profit_rates = 0:step:round(R, digits = precision)

    env = LPEnvelope(A, profit_rates, effects_sectors)

    # FIXME: Das gehört hier nicht her
    # Set l
    map((x, y) -> set_objective_coefficient(model_intensities, x, y),
        model_intensities[:x], l)
    # This is necessary because of usa1982
    if !all(d .>= 0)
        d = ones(env.n_goods)
    end
    # Set d
    set_normalized_rhs.(model_intensities[:d], d)

    while R < R_max
        R_max = -1.0
        for (i, r) in enumerate(profit_rates)
            verbose && println("Profit rate: $r")
            # Get the indices for the columns with the highest intensity per sector
            function highest_intensity_indices(sector)
                max_index = argmax(j -> env.intensities[j, i], sector:env.n_goods:length(l))
                return max_index
            end
    
            # von Neumann Ansatz
            C = (B - (1 + r) * A)
            # See Han p. 169
            qq = ones(size(A, 2))
            if all(C * qq .- d .< 0.001)
                println("B-(1+$r) equals final demand d for r=$r")
            end
    
            # Compute intensities
            modify_r!(model_intensities, r)
            optimize!(model_intensities) # min l * x s.t. C * x ≥ d AND x .≥ 0
            error_msg = "$(termination_status(model_intensities))\nProbably unfeasible d!"
            @assert is_solved_and_feasible(model_intensities) error_msg
            env.intensities[:, i] = replace_with_zero.(value.(model_intensities[:x]))
            env.chosen_technology[:, i] = [highest_intensity_indices(sector)
                                      for sector in 1:env.n_goods]
        end

        # We compute the truncated stuff and prices only at switch points
        for (i, tech) in enumerate(eachcol(env.chosen_technology)[begin:(end - 1)])
            # Check whether we are at switch point
            switch = tech[1:effects_sectors] != env.chosen_technology[1:effects_sectors, i + 1]
            if switch
                for j in i:(i + 1)
                    # Compute truncated intensities
                    d = ones(env.n_goods)
                    A_trunc = view(A, :, env.chosen_technology[:, j])  # filter columns
                    l_trunc = vec(view(l, :, env.chosen_technology[:, j]))'  # filter columns
                    B_trunc = I(size(A_trunc, 1))  # initialize identity matrix
                    C_trunc = B_trunc - A_trunc
                    set_normalized_rhs.(model_intensities_trunc[:con], d)
                    modify_C!(model_intensities_trunc, model_intensities_trunc[:x], C_trunc)
                    optimize!(model_intensities_trunc)  # min l_trunc * x s.t. C_trunc * x ≥ d AND x .≥ 0
                    @assert is_solved_and_feasible(model_intensities_trunc)
                    env.intensities_trunc[:, j] = value.(model_intensities_trunc[:x])

                    # Compute prices
                    r = profit_rates[j]
                    C = (B - (1 + r) * A)
                    C_trunc = view(C, :, env.chosen_technology[:, j])  # filter columns
                    set_normalized_rhs.(model_prices[:con], vec(l_trunc))
                    set_objective_coefficient.(model_prices, model_prices[:p], d)
                    modify_C!(model_prices, model_prices[:p], permutedims(C_trunc))
                    optimize!(model_prices)  # max p * d s.t. C' * p ≤ l AND p .≥ 0
                    @assert is_solved_and_feasible(model_prices)
                    env.prices[:, j] = value.(model_prices[:p])

                    # Compute (l * x)
                    env.lx[j] = l_trunc[:, 1:effects_sectors] * env.intensities_trunc[:, j][1:effects_sectors]
                    # Compute (pA)'
                    env.pA[:, j] = vec(A_trunc[1:effects_sectors, 1:effects_sectors]' * env.prices[:, i][1:effects_sectors])
                    # Compute A * x / (l * x)
                    env.right_side_factor[:, j] = A_trunc[1:effects_sectors, 1:effects_sectors] *
                                                env.intensities_trunc[:, j][1:effects_sectors] / env.lx[j]
                end

                # Compute p * A * x / (l * x)
                push!(env.capital_intensities,
                    [
                        profit_rates[i] => env.prices[:, i][1:effects_sectors]' * env.right_side_factor[:, i],
                        profit_rates[i + 1] => env.prices[:, i][1:effects_sectors]' * env.right_side_factor[:, i + 1]
                    ]
                )
                push!(env.lx_at_switch, [
                        profit_rates[i] => env.lx[i],
                        profit_rates[i + 1] => env.lx[i + 1]
                    ]
                )
                push!(env.pA_at_switch,
                    [
                        profit_rates[i] => env.pA[:, i][1:effects_sectors],
                        profit_rates[i + 1] => env.pA[:, i + 1][1:effects_sectors]
                    ]
                )
                push!(env.l_at_switch,
                    [
                        profit_rates[i] => vec(view(l, :, env.chosen_technology[:, i]))[1:effects_sectors],
                        profit_rates[i + 1] => vec(view(l, :, env.chosen_technology[:, i + 1]))[1:effects_sectors]
                    ]
                )
                push!(env.intensities_at_switch,
                [
                    profit_rates[i] => env.intensities_trunc[:, i][1:effects_sectors],
                    profit_rates[i + 1] => env.intensities_trunc[:, i + 1][1:effects_sectors]
                ]
                )
                push!(env.prices_switch, profit_rates[i] => env.prices[:, i][1:effects_sectors])
                push!(env.technologies_switch,
                    [profit_rates[i] => env.chosen_technology[:, i][1:effects_sectors],
                    profit_rates[i + 1] => env.chosen_technology[:, i + 1][1:effects_sectors]]
                )
            end
        end
    end

    profit_rates_to_names = Dict(profit_rates .=>
        string.(round.(profit_rates, digits = precision)))
    df_intensities = DataFrame(
        env.intensities, map(r -> profit_rates_to_names[r], profit_rates))

    if save_all
        switches = Dict(
            "capital_intensities" => env.capital_intensities,
            "pA" => env.pA_at_switch,
            "lx" => env.lx_at_switch,
            "l" => env.l_at_switch,
            "prices" => env.prices_switch,
            "intensities" => env.intensities_at_switch,
            "technology" => env.technologies_switch
        )
    else
        switches = Dict(
            "capital_intensities" => env.capital_intensities,
            "l" => env.l_at_switch,
            "prices" => env.prices_switch,
            "intensities" => env.intensities_at_switch,
            "technology" => env.technologies_switch
        )
    end

    return df_intensities, profit_rates_to_names, profit_rates, switches
end