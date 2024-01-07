module GaussianOptics

export
    # Types
    GaussianBeamParams,
    GaussianBeamAngle, 
    GaussianBeamWaist, 
    # Functions
    beamWaist,
    divergentAngle,
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

#=
                     ┌─────────────────┐
                     │ Rayleigh Length │
                     │        zᵣ       │
                     └─────────────────┘
                     ╱                 ╲
        zᵣ = πw₀²/λ ╱                   ╲ zᵣ = 4λ/πΘ²
                   ╱                     ╲
  ┌────────────────┐      Θ = 2λ/πw₀     ┌─────────────────┐
  │ Min. Beamwaist │─────────────────────│ Divergent Angle │
  │       w₀       │─────────────────────│        Θ        │
  └────────────────┘     w₀ = 2λ/πΘ      └─────────────────┘
=#

beamWaist(w::GaussianBeamWaist)               = w.waist
beamWaist(λ::Real, Θ::GaussianBeamAngle)      = 0.6366197723675814 * λ / Θ.angle

divergentAngle(Θ::GaussianBeamAngle)          = Θ.angle
divergentAngle(λ::Real, w::GaussianBeamWaist) = 0.6366197723675814 * λ / w.waist

rayleighLength(λ::Real, Θ::GaussianBeamAngle) = 1.2732395447351628 * λ / abs2(Θ.angle)
rayleighLength(λ::Real, w::GaussianBeamWaist) = π * abs2(w.waist) / λ

end # module GaussianOptics
