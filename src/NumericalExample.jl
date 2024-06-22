module NumericalExample

export model_parameters, calculate_steady_state, tax_revenue, period_budget, euler_residual,
    labor_FOC_residual, approximation_setup, calculate_initial_guess, approximated_functions,
    residuals_callable, solve_model, default_colorscheme, plot_vs_time

using ArgCheck: @argcheck
import ColorSchemes
using DocStringExtensions: FIELDS, SIGNATURES
using LogExpFunctions: logistic
using InverseFunctions: inverse
using Miter
using SpectralKit: Chebyshev, InteriorGrid, SemiInfRational, dimension, grid,
    linear_combination
using TrustRegionMethods: ForwardDiff_wrapper, trust_region_solver, SolverStoppingCriterion

####
#### model and solution conditions
####

"""
A container for the model parameters. Use [`model_parameters`](@ref) to create, which
performs some simple checks.

$(FIELDS)
"""
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

function Base.show(io::IO, model::T) where {T<:ModelParameters}
    print(io, "NumericalExample.model_parameters(; ")
    join(io, (string(k) * " = " * repr(getfield(model, k)) for k in fieldnames(T)), ", ")
    print(io, ")")
end

"""
$(SIGNATURES)

Create a value that holds the parameters of the model.

# Keyword arguments

- `θ`: the parameter ``θ \\log(c) + (1-θ) \\log(1-ℓ)`` in the utility function
- `β`: the discount factor
- `α`: the capital share, in the production function ``k^α ℓ^{1-α}``
- `δ`: the depreciation rate
- `τ_k`, `τ_ℓ`: tax rates for capital and labor, respectively
"""
function model_parameters(; θ = 0.5, β = 0.95, α = 1/3, δ = 0.05, τ_k = 0.2, τ_ℓ = 0.2)
    @argcheck 0 < θ < 1
    @argcheck 0 < β < 1
    @argcheck 0 < α < 0.98 "Invalid α. Current solver does not handle α ≈ 0.99 nicely, disabled until implemented."
    @argcheck 0 < δ < 1
    @argcheck 0 < τ_k < 1
    @argcheck 0 < τ_ℓ < 1
    ModelParameters(; θ, β, α, δ, τ_k, τ_ℓ)
end

"""
$(SIGNATURES)

Calculate the rent and the wage in the `model`, using capital `k` and labor `ℓ`.
"""
function rent_and_wage(model::ModelParameters, k, ℓ)
    (; α) = model
    @argcheck k > 0
    @argcheck 0 < ℓ ≤ 1
    ratio = k / ℓ
    r = α * ratio^(α - 1)
    w = (1 - α) * ratio^α
    r, w
end

"""
$(SIGNATURES)

Calculate the tax revenue in the `model`, using capital `k` and labor `ℓ`.
"""
function tax_revenue(model::ModelParameters, k, ℓ)
    (; τ_k, τ_ℓ) = model
    r, w = rent_and_wage(model, k, ℓ)
    τ_k * r * k + τ_ℓ * r * ℓ
end

"""
$(SIGNATURES)

Calculate the budget of the agent each period (can be spent on consumption and
next-period capital).
"""
function period_budget(model::ModelParameters, k, ℓ)
    (; τ_k, τ_ℓ, δ) = model
    r, w = rent_and_wage(model, k, ℓ)
    ((1 - τ_k) * r + (1 - δ)) * k + (1 - τ_ℓ) * w * ℓ
end

"""
$(SIGNATURES)

Calculate the residual of the intertemporal equilibrium condition (Euler equation) at a
given point in time.

`c` is consumption for this period, `k′`, `ℓ′`, `c′` are capital, labor, and consumption
*next* period, respectively.
"""
function euler_residual(model::ModelParameters; c, k′, ℓ′, c′)
    (; β, θ, τ_k, δ) = model
    r′, _ = rent_and_wage(model, k′, ℓ′)
    lhs = θ / c
    rhs = β * θ / c′ * ((1 - τ_k) * r′ + 1 - δ)
    rhs / lhs - 1
end

"""
$(SIGNATURES)

Calculate the residual of the intratemporal first order condition.

 `k`, `ℓ`, `c` are capital, labor, and consumption this period, respectively.
"""
function labor_FOC_residual(model::ModelParameters; k, ℓ, c)
    (; θ, τ_ℓ) = model
    _, w = rent_and_wage(model, k, ℓ)
    θ / c * (1 - τ_ℓ) * w - (1-θ) / (1-ℓ)
end

"""
$(SIGNATURES)

Calculate the steady state solution, returning a container with fields `c̄`, `k̄`, `ℓ̄`.
"""
function calculate_steady_state(model::ModelParameters)
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

####
#### approximate solution calculation
####

"""
A container for the approximation setup. Use [`approximaton_setup`](@ref) to create.

$(FIELDS)
"""
struct ApproximationSetup{B,TP,TU}
    basis::B
    positive_transformation::TP
    unit_transformation::TU
end

"""
$(SIGNATURES)

Return a value that contains all parameters for the approximation.

- `basis0`: the untransformed basis (with a finite domain), defaults to a Chebyshev
  approximation with 10 elements, nodes on the interior of the domain

- `mid_t`: the time period for more or less the middle node. Approximately, it should be
  the time horizon we are interested in, but the method is fairly accurate without
  determining this precisely

- `positive_transformation`, `unit_transformation`: functions that transform ``ℝ`` to
  ``(0, ∞)`` and ``(0, 1)``. The defaults are `exp` and `logistic`, and you probably do
  not need to change these.
"""
function approximation_setup(;
                             basis0 = Chebyshev(InteriorGrid(), 10),
                             mid_t = 5.0,
                             positive_transformation = exp,
                             unit_transformation = logistic)
    basis = basis0 ∘ SemiInfRational(0, mid_t)
    ApproximationSetup(basis, positive_transformation, unit_transformation)
end

function calculate_initial_guess(approx::ApproximationSetup, steady_state)
    (; c̄, k̄, ℓ̄) = steady_state
    (; basis, positive_transformation, unit_transformation) = approx
    D = dimension(approx.basis)
    θ = zeros(3 * D)
    Ip = inverse(positive_transformation)
    θ[1] = Ip(k̄)
    θ[D + 1] = Ip(c̄)
    θ[2 * D + 1] = inverse(unit_transformation)(ℓ̄)
    θ
end

function approximated_functions(approx::ApproximationSetup, θ)
    (; basis, positive_transformation, unit_transformation) = approx
    D = dimension(basis)
    @argcheck length(θ) == D * 3
    θ3 = reshape(θ, D, 3)
    k̃ = positive_transformation ∘ linear_combination(basis, θ3[:, 1])
    c̃ = positive_transformation ∘ linear_combination(basis, θ3[:, 2])
    ℓ̃ = unit_transformation ∘ linear_combination(basis, θ3[:, 3])
    (; k̃, c̃, ℓ̃)
end

function calculate_residuals(model::ModelParameters, approx::ApproximationSetup,
                             steady_state, k0, approxf)
    (; c̃, k̃, ℓ̃) = approxf
    (; c̄, k̄, ℓ̄) = steady_state
    residuals = [c̃(Inf) / c̄ - 1, k̃(Inf) / k̄ - 1, ℓ̃(Inf) / ℓ̄ - 1, k̃(0) / k0 - 1]
    for t in grid(approx.basis)
        c = c̃(t)
        k = k̃(t)
        ℓ = ℓ̃(t)
        c′ = c̃(t + 1)
        k′ = k̃(t + 1)
        ℓ′ = ℓ̃(t + 1)
        push!(residuals,
              period_budget(model, k, ℓ) / (c + k′) - 1,
              euler_residual(model; c, k′, ℓ′, c′),
              labor_FOC_residual(model; c, ℓ, k))
    end
    residuals
end

Base.@kwdef struct ResidualsCallable{M,A,S,T<:Real}
    model::M
    approx::A
    steady_state::S
    k0::T
end

function residuals_callable(model, approx, k0)
    steady_state = calculate_steady_state(model)
    ResidualsCallable(; model, approx, steady_state, k0)
end

function calculate_initial_guess(rc::ResidualsCallable)
    calculate_initial_guess(rc.approx, rc.steady_state)
end

function (rc::ResidualsCallable)(θ)
    (; model, approx, steady_state, k0) = rc
    approxf = approximated_functions(approx, θ)
    calculate_residuals(model, approx, steady_state, k0, approxf)
end

Base.@kwdef struct ModelSolution{TC,TK,TL}
    converged::Bool
    c̃::TC
    k̃::TK
    ℓ̃::TL
end

function Base.show(io::IO, solution::ModelSolution)
    print(io, "model solution ")
    if solution.converged
        printstyled(io, "converged"; color = :green)
    else
        printstyled(io, "did not converge"; color = :red)
    end
end

function solve_model(model::ModelParameters, k0;
                     approx = approximation_setup(),
                     tol = 1e-2)
    rc = residuals_callable(model, approx, k0)
    θ0 = calculate_initial_guess(rc)
    sol = trust_region_solver(ForwardDiff_wrapper(rc, length(θ0)), θ0;
                              stopping_criterion = SolverStoppingCriterion(tol))
    (; converged, x) = sol
    approxf = approximated_functions(approx, x)
    ModelSolution(; converged, approxf...)
end

####
#### plotting code
####

default_colorscheme() = ColorSchemes.Dark2_8

function plot_vs_time(fs, ts, label = ""; colors = default_colorscheme(),
                      graph_labels = nothing, label_relative_offset = 0.05)
    values = mapreduce(fs, hcat, ts)
    label_offset = -(extrema(values)...) * label_relative_offset
    function _graph(v, color, i)
        p = Any[Lines(zip(ts, v); color)]
        if graph_labels ≠ nothing
            j = div(length(v), 2)
            push!(p, Annotation((ts[j], v[j] - label_offset), textcolor(color, graph_labels[i]);
                                bottom = true))
        end
        p
    end
    graphs = mapreduce(_graph, vcat, eachrow(values), colors, axes(values, 1))
    Plot(graphs,
         x_axis = Axis.Linear(; label = "time"), y_axis = Axis.Linear(; label))
end

end # module
