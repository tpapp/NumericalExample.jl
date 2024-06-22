# Numerical approximation

## Defining the model

We set up a container for the model parameters, but it is not exported from the module.
```@docs
NumericalExample.ModelParameters
```
We use a function to instantiate it, which is better style in Julia and allows for seamless extension of the interface later.
```@docs
model_parameters
```

## Coding the relevant equilibrium conditions

We code the elements of the model in simple functions, which are easier to check and debug.
```@docs
NumericalExample.rent_and_wage
NumericalExample.tax_revenue
NumericalExample.period_budget
```

We encapsulate equilibrium conditions into two simple functions.
```@docs
NumericalExample.euler_residual
NumericalExample.labor_FOC_residual
```

## Steady state calculation

Then we code the steady state as derived in the previous section.
```@docs
calculate_steady_state
```

## Spectral approximation

We solve the full model in the time domain, using Chebyshev polynomials transformed to ``[0, \infty]`` [boyd2001chebyshev](@cite). 

The method is similar to [judd_parametric_2002](@cite), except that we find that rational transformations of Chebyshev polynomials work well without building in exponential decay. Also, in contrast to a lot of texts, we approximate not only ``k(t)`` but also ``c(t)`` and ``\ell(t)``, instead of deriving them from equilibrium conditions. This leads to a larger parameter space for the approximation, but ensures that the residuals are always well-defined and have finite derivatives, leading the solver quickly to the equilibrium. This requires that we transform values to be positive (for ``c(t)`` and ``k(t)``) and to be within ``(0, 1)`` (for ``\ell(t)``).

We first define the approximation setup.

```@docs
approximation_setup
```

Spectral approximations usually work best if we have a reasonably good initial guess. For this we use the steady state solution ``\forall t > 0``, since that is always a solution with ``k(0) = k_0`` for any initial capital level ``k_0``.

```@docs
calculate_initial_guess
```

## Solving the model

To solve the model, we define a mapping from the approximation parameters to the residuals, and use a trust region method to solve for the roots, with automatic differentiation.

The solver is encapsulated into the function below.

```@docs
solve_model
```
