module GaussianOptics

export
    # Types
    GaussianBeamParams,
    GaussianBeamAngle, 
    GaussianBeamWaist, 
    # Functions
    rayleighLength,
    # Constants
    GBParams,
    GBAngle,
    GBWaist

abstract type GaussianBeamParams end

const GBParams = GaussianBeamParams

struct GaussianBeamAngle <: GaussianBeamParams angle::Float64 end
struct GaussianBeamWaist <: GaussianBeamParams waist::Float64 end

const GBAngle = GaussianBeamAngle
const GBWaist = GaussianBeamWaist

rayleighLength(λ::Real, Θ::GaussianBeamAngle) = 1.2732395447351628 * λ / abs2(Θ.angle)
rayleighLength(λ::Real, w::GaussianBeamWaist) = π * abs2(w.waist) / λ

end # module GaussianOptics
