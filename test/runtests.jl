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
        (; k̄, ℓ̄, c̄) = steady_state(model)
        fail = false
        fail |= (@test period_budget(model, k̄, ℓ̄) ≈ k̄ + c̄ atol = 1e-8) == Test.Fail
        @test euler_residual(model; c = c̄, c′ = c̄, k′ = k̄, ℓ′ = ℓ̄) ≈ 0 atol = 1e-8
        @test labor_FOC_residual(model; c = c̄, k = k̄, ℓ = ℓ̄) ≈ 0 atol = 1e-8
        if fail
            @show model
        end
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
