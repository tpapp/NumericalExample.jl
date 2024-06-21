using NumericalExample
using Test

function random_parameters()
    θ = rand()
    β = rand() * 0.2 + 0.8
    α = rand() * 0.98
    δ = rand()
    τ_k = rand()
    τ_ℓ = rand()
    model_parameters(; θ, β, α, δ, τ_k, τ_ℓ)
end

@testset "test steady state solver" begin
    for _ in 1:100
        model = random_parameters()
        (; k̄, ℓ̄, c̄) = calculate_steady_state(model)
        @test period_budget(model, k̄, ℓ̄) ≈ k̄ + c̄ atol = 1e-8 rtol = 1e-8
        @test euler_residual(model; c = c̄, c′ = c̄, k′ = k̄, ℓ′ = ℓ̄) ≈ 0 atol = 1e-8
        @test labor_FOC_residual(model; c = c̄, k = k̄, ℓ = ℓ̄) ≈ 0 atol = 1e-8
    end
end

@testset "initial guess calculations" begin
    approx = approximation_setup()
    for _ in 1:100
        model = random_parameters()
        steady_state = calculate_steady_state(model)
        s0 = calculate_initial_guess(approx, steady_state)
        f = approximated_functions(approx, s0)
        for t in vcat(range(0.1, 10.0, 50), [0, Inf])
            @test f.k̃(t) ≈ steady_state.k̄
            @test f.c̃(t) ≈ steady_state.c̄
            @test f.ℓ̃(t) ≈ steady_state.ℓ̄
        end
    end
end

@testset "residuals callable calculations with steady state" begin
    approx = approximation_setup()
    for _ in 1:100
        model = random_parameters()
        steady_state = calculate_steady_state(model)
        rc = residuals_callable(model, approx, steady_state.k̄)
        θ0 = @inferred calculate_initial_guess(rc)
        r0 = @inferred rc(θ0)
        @test sum(abs2, r0) ≤ 1e-8
    end
end

@testset "solving for random parameters and steady states" begin
    for _ in 1:50
        model = random_parameters()
        sol = solve_model(model, calculate_steady_state(model).k̄ * (rand() + 0.5))
        @test sol.converged
    end
end

## NOTE add JET to the test environment, then uncomment
# using JET
# @testset "static analysis with JET.jl" begin
#     @test isempty(JET.get_reports(report_package(NumericalExample, target_modules=(NumericalExample,))))
# end

## NOTE add Aqua to the test environment, then uncomment
# @testset "QA with Aqua" begin
#     import Aqua
#     Aqua.test_all(NumericalExample; ambiguities = false)
#     # testing separately, cf https://github.com/JuliaTesting/Aqua.jl/issues/77
#     Aqua.test_ambiguities(NumericalExample)
# end
