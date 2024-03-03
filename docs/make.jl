using Documenter, Documenter.Remotes
using ITensors, PseudomodesTTEDOPA

makedocs(
    sitename = "PseudomodesTTEDOPA",
    format = Documenter.HTML(),
    repo = Remotes.GitHub("phaerrax", "PseudomodesTTEDOPA.jl"),
    modules = [PseudomodesTTEDOPA],
    checkdocs = :none,
    draft = true,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
