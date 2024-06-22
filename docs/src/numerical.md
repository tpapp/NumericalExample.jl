# Numerical approximation

We set up a container for the model parameters, but it is not exported from the module.
```@docs
NumericalExample.ModelParameters
```
We use a function to instantiate it, which is better style in Julia and allows for seamless extension of the interface later.
```@docs
model_parameters
```

We code the elements of the model in simple functions, which are easier to check and debug.
```@docs
NumericalExample.rent_and_wage
NumericalExample.tax_revenue
NumericalExample.period_budget
```

```@docs
NumericalExample.euler_residual
NumericalExample.labor_FOC_residual
```
