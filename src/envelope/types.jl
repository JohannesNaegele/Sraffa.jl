Base.@kwdef struct LPEnvelope{S, T, U, EAM, EAI, EAV, VP, VPV, VPI, P} <: EnvelopeMethod
    A::S
    l::T
    d::U
    # Number of goods
    n_goods::Int64
    n_countries::Int64

    # Activity levels (q)
    intensities::EAM
    # Truncated activity levels; calculated via C = (B - A)
    intensities_trunc::EAM
    pA::EAM
    right_side_factor::EAM
    lx::EAV
    # Prices computed from truncated intensities via the dual representation
    prices::EAM
    # Indices of chosen technologies; chosen_technology[i] ∈ i * n_goods
    chosen_technology::EAI

    # These are data structures used to store only at switches
    # Capital intensities (pAx / lx)
    capital_intensities::VP
    lx_at_switch::VP
    intensities_at_switch::VPV
    # Prices (specified to be exactly the same before/after switch; it doesn't make a difference though)
    prices_switch::P
    technologies_switch::VPI
    pA_at_switch::VPV
    l_at_switch::VPV
end


function LPEnvelope(A, l, d, profit_rates, effects_sectors)
    n_goods = size(A, 1)
    LPEnvelope(
        A = A,
        l = l,
        d = d,
        n_goods = n_goods,
        n_countries = div(size(A, 2), n_goods),
        intensities = ElasticArray{Float64}(undef, size(A, 2), length(profit_rates)),
        intensities_trunc = ElasticArray{Float64}(undef, n_goods, length(profit_rates)),
        pA = ElasticArray{Float64}(undef, length(effects_sectors), length(profit_rates)),
        right_side_factor = ElasticArray{Float64}(undef, length(effects_sectors), length(profit_rates)),
        lx = ElasticArray{Float64}(undef, length(profit_rates)),
        prices = ElasticArray{Float64}(undef, n_goods, length(profit_rates)),
        chosen_technology = ElasticArray{Int64}(undef, n_goods, length(profit_rates)),
        capital_intensities = Vector{Pair{Float64, Float64}}[],
        lx_at_switch = Vector{Pair{Float64, Float64}}[],
        intensities_at_switch = Vector{Pair{Float64, Vector{Float64}}}[],
        prices_switch = Pair{Float64, Vector{Float64}}[],
        technologies_switch = Vector{Pair{Float64, Vector{Int64}}}[],
        pA_at_switch = Vector{Pair{Float64, Vector{Float64}}}[],
        l_at_switch = Vector{Pair{Float64, Vector{Float64}}}[]
    )
end

function extend!(env::LPEnvelope, steps)
    append!(env.intensities, zeros(size(env.intensities, 1), steps))
    append!(env.intensities_trunc, zeros(size(env.intensities_trunc, 1), steps))
    append!(env.pA, zeros(size(env.pA, 1), steps))
    append!(env.right_side_factor, zeros(size(env.right_side_factor, 1), steps))
    append!(env.lx, zeros(steps))
    append!(env.prices, zeros(size(env.prices, 1), steps))
    append!(env.chosen_technology, zeros(size(env.chosen_technology, 1), steps))
end

function save(env::LPEnvelope, save_all)
    if save_all
        return Dict(
            "capital_intensities" => env.capital_intensities,
            "pA" => env.pA_at_switch,
            "lx" => env.lx_at_switch,
            "l" => env.l_at_switch,
            "prices" => env.prices_switch,
            "intensities" => env.intensities_at_switch,
            "technology" => env.technologies_switch
        )
    else
        return Dict(
            "capital_intensities" => env.capital_intensities,
            "l" => env.l_at_switch,
            "prices" => env.prices_switch,
            "intensities" => env.intensities_at_switch,
            "technology" => env.technologies_switch
        )
    end
end
