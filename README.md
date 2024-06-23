# NumericalExample.jl

![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
[![build](https://github.com/tpapp/NumericalExample.jl/workflows/CI/badge.svg)](https://github.com/tpapp/NumericalExample.jl/actions?query=workflow%3ACI)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://tpapp.github.io/NumericalExample.jl/dev)
<!-- [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) -->

## Overview

This a demonstration of numerical model solving and a basic policy analysis environment in Julia.

The emphasis is on showcasing good practices for

1. software documentation, including docstrings in the [source code](./src/NumericalExample.jl) and the [generated documentation](https://tpapp.github.io/NumericalExample.jl/dev)
2. unit [tests](./test/runtest.jl) that ensure the code is always in a correct state,
3. some interactive analysis and plotting,
4. scientific software development.

Because of the focus on the above, the model is kept trivially simple. An actual economic model would be more complex and be calibrated/estimated using data (which here is omitted for simplicity), but code would be organized along the same principles.

## Installation 

This is not registered package. Install in an interactive environment with

```julia
julia> ]add https://github.com/tpapp/NumericalExample.jl
```

Julia's package delivery mechanism ensures that all required packages are installed. The [manifest](./Manifest.toml) is commiteed to the repository so you get the exact same versions --- this is important for scientific code.
