module GaussianOptics

export
    # Types
    GaussianBeamParams,
    GaussianBeamAngle, 
    GaussianBeamWaist, 
    GaussianBeamState,
    TransferOperator,
    ConstSpace,
    ThinLens,
    # Functions
    beamWaist,
    divergentAngle,
    rayleighLength,
    # Constants
    GBParams,
    GBAngle,
    GBWaist,
    GBState

abstract type GaussianBeamParams end

const GBParams = GaussianBeamParams

struct GaussianBeamAngle <: GaussianBeamParams angle::Float64    end
struct GaussianBeamWaist <: GaussianBeamParams waist::Float64    end
struct GaussianBeamState <: GaussianBeamParams state::ComplexF64 end

GaussianBeamState(λ::Real, p::GaussianBeamParams; z::Real=0.0) = GaussianBeamState(complex(z, -rayleighLength(λ, p)))

const GBAngle = GaussianBeamAngle
const GBWaist = GaussianBeamWaist
const GBState = GaussianBeamState

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
beamWaist(λ::Real, p::GaussianBeamState)      = sqrt(λ / (π * imag(inv(p.state))))

divergentAngle(Θ::GaussianBeamAngle)          = Θ.angle
divergentAngle(λ::Real, w::GaussianBeamWaist) = 0.6366197723675814 * λ / w.waist

rayleighLength(λ::Real, Θ::GaussianBeamAngle) = 1.2732395447351628 * λ / abs2(Θ.angle)
rayleighLength(λ::Real, w::GaussianBeamWaist) = π * abs2(w.waist) / λ

# Ray Transfer Matrix Analysis

abstract type TransferOperator end

struct ConstSpace <: TransferOperator d::Float64 end
struct ThinLens   <: TransferOperator f::Float64 end

constSpace(q::ComplexF64, d::Real) = q + d
thinLens(q::ComplexF64, f::Real)   = q / (1.0 - q / f)

Base.:*(op::ConstSpace, s::GaussianBeamState) = GaussianBeamState(constSpace(s.state, op.d))
Base.:*(op::ThinLens,   s::GaussianBeamState) = GaussianBeamState(  thinLens(s.state, op.f))

Base.:*(s::GaussianBeamState, op::ConstSpace) = GaussianBeamState(constSpace(s.state, op.d))
Base.:*(s::GaussianBeamState, op::ThinLens)   = GaussianBeamState(  thinLens(s.state, op.f))

end # module GaussianOptics
