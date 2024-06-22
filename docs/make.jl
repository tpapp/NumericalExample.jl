# see documentation at https://juliadocs.github.io/Documenter.jl/stable/

using Documenter, DocumenterCitations, NumericalExample

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

makedocs(;
         modules = [NumericalExample],
         format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true",
                                  assets = ["assets/custom.css"]),
         authors = "Tamás K. Papp",
         sitename = "NumericalExample.jl",
         pages = Any["index.md",
                     "theory.md",
                     "numerical.md",
                     "application.md"],
         plugins = [bib],
         clean = true,
         checkdocs = :exports,
         warnonly = true # Documenter.except(:missing_docs),
         )

# Some setup is needed for documentation deployment, see “Hosting Documentation” and
# deploydocs() in the Documenter manual for more information.
deploydocs(
    repo = "github.com/tpapp/NumericalExample.jl.git",
    push_preview = true
)
