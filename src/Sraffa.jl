module Sraffa

export compute_envelope

# Helper functions
include("helpers.jl")

# Envelope calculation
include("envelope/lp.jl")
include("envelope/binary.jl")
include("envelope/vfz.jl")
include("envelope/envelope.jl")

end