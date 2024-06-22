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

## Numerical approximation

We solve the full model in the time domain, using Chebyshev polynomials transformed to ``[0, \infty]`` [boyd2001chebyshev](@cite). The method is similar to [judd_parametric_2002](@cite), except that we find that rational transformations of Chebyshev polynomials work well without building in exponential decay.

We first define the approximation setup.

```@docs
approximation_setup
```
