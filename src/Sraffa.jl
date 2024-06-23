module Sraffa

export compute_envelope
export create_intensities_r_solver, create_intensities_solver, create_prices_solver
export modify_A!, modify_C!, modify_r!
export real_eigvals, compute_R, replace_with_zero, safe_divide
export n_switches, n_reswitches, switch_cases, switch_info
export compute_w, plot_wage_curves

using Logging
using LinearAlgebra
using DataFrames
using ElasticArrays
using LazyArrays
using StrideArrays
using JuMP
using RecipesBase
using FLoops
using LoopVectorization
using FastBroadcast

# Helper functions
include("helpers.jl")

# Envelope calculation
abstract type EnvelopeMethod end

include("envelope/types.jl")
include("envelope/lp.jl")
include("envelope/piecewise.jl")
include("envelope/vfz.jl")
include("envelope/envelope.jl")

# Evaluation of results
include("evaluate.jl")

end # module Sraffa