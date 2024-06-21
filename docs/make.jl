# see documentation at https://juliadocs.github.io/Documenter.jl/stable/

using Documenter, NumericalExample

# mathengine = Documenter.MathJax()

makedocs(
    modules = [NumericalExample],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true",
                             assets = ["assets/custom.css"], #=, mathengine=#),
    authors = "Tamás K. Papp",
    sitename = "NumericalExample.jl",
    pages = Any["index.md",
                "theory.md",
                "numerical.md",
                "application.md"],
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

# Some setup is needed for documentation deployment, see “Hosting Documentation” and
# deploydocs() in the Documenter manual for more information.
deploydocs(
    repo = "github.com/tpapp/NumericalExample.jl.git",
    push_preview = true
)
