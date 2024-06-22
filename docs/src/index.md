# Overview

This is a self-contained worked example of solving a *very simple* macroeconomic model in Julia and doing some comparative analysis with the results.

Since this is just a demonstration, the model is of course not realistic, more like something you would see in a homework problem in an introductory graduate course. The model would need to be extended with a more detailed production structure, capital and labor frictions, nontrivial consumers, etc, and estimated based on data, to be considered useful for actual policy analysis. For an introductory text, see [krusell2014real](@cite) and [acemoglu2008introduction](@cite).

The [code](https://github.com/tpapp/NumericalExample.jl/blob/master/src/NumericalExample.jl) is kept simple and it is well-documented. An important part of this demonstration is how all components are [tested](https://github.com/tpapp/NumericalExample.jl/blob/master/test/runtests.jl) in detail to ensure correctness.

## Bibliography

```@bibliography
```
