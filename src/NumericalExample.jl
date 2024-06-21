"""
Placeholder for a short summary about NumericalExample.
"""
module NumericalExample

export model_parameters, steady_state, period_budget, euler_residual, labor_FOC_residual

using ArgCheck: @argcheck

Base.@kwdef struct ModelParameters{T}
    "weight on (log) consumption in period utility"
    θ::T
    "discount factor"
    β::T
    "capital share"
    α::T
    "depreciation"
    δ::T
    "tax rate on capital returns"
    τ_k::T
    "tax rate on wages"
    τ_ℓ::T
end

function model_parameters(; θ = 0.5, β = 0.95, α = 1/3, δ = 0.05, τ_k = 0.2, τ_ℓ = 0.2)
    @argcheck 0 < θ < 1
    @argcheck 0 < β < 1
    @argcheck 0 < α < 0.98 "Invalid α. Current solver does not handle α ≈ 0.99 nicely, disabled until implemented."
    @argcheck 0 < δ < 1
    @argcheck 0 < τ_k < 1
    @argcheck 0 < τ_ℓ < 1
    ModelParameters(; θ, β, α, δ, τ_k, τ_ℓ)
end

function rent_and_wage(model::ModelParameters, k, ℓ)
    (; α) = model
    @argcheck k > 0
    @argcheck 0 < ℓ ≤ 1
    ratio = k / ℓ
    r = α * ratio^(α - 1)
    w = (1 - α) * ratio^α
    r, w
end

function period_budget(model::ModelParameters, k, ℓ)
    (; τ_k, τ_ℓ, δ) = model
    r, w = rent_and_wage(model, k, ℓ)
    ((1 - τ_k) * r + (1 - δ)) * k + (1 - τ_ℓ) * w * ℓ
end

function euler_residual(model::ModelParameters; c, k′, ℓ′, c′)
    (; β, θ, τ_k, δ) = model
    r′, _ = rent_and_wage(model, k′, ℓ′)
    lhs = θ / c
    rhs = β * θ / c′ * ((1 - τ_k) * r′ + 1 - δ)
    rhs / lhs - 1
end

function labor_FOC_residual(model::ModelParameters; k, ℓ, c)
    (; θ, τ_ℓ) = model
    _, w = rent_and_wage(model, k, ℓ)
    θ/c*(1 - τ_ℓ)*w - (1-θ)/(1-ℓ)
end

function steady_state(model::ModelParameters)
    (; α, β, θ, δ, τ_k, τ_ℓ) = model
    λ1mα = (1/β - (1 - δ)) / ((1 - τ_k) * α) # λ^{1-α}
    c̄_k̄_ratio = (1 - τ_k) * α * λ1mα + (1 - δ) + (1- τ_ℓ) * (1 - α) * λ1mα - 1
    λ = λ1mα^(1 / (1-α))
    R = θ / (c̄_k̄_ratio) * (1- τ_ℓ) * (1-α) * λ^(-α) # RHS of steady state intratemporal
    k̄ = 1 / ((1 - θ) / R + λ)
    ℓ̄ = λ * k̄
    c̄ = c̄_k̄_ratio * k̄
    (; c̄, k̄, ℓ̄)
end

end # module
