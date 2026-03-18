using Documenter
using FrontTrackingMethods

makedocs(
    modules = [FrontTrackingMethods],
    authors = "PenguinxCutCell contributors",
    sitename = "FrontTrackingMethods.jl",
    format = Documenter.HTML(
        canonical = "https://PenguinxCutCell.github.io/FrontTrackingMethods.jl",
        repolink = "https://github.com/PenguinxCutCell/FrontTrackingMethods.jl",
        collapselevel = 2,
    ),
    pages = [
        "Home" => "index.md",
        "Benchmarks" => "benchmarks.md",
        "Diagnostics" => "diagnostics.md",
        "Remeshing" => "remeshing.md",
        "Topology change (v0.3)" => "topology_change.md",
        "API Reference" => "api.md",
        "Examples" => "examples.md",
        "Plotting and animations" => "plotting.md",
        "Installation" => "installation.md",
        "Usage" => "usage.md",
        "Testing" => "testing.md"
        ],
    pagesonly = true,
    warnonly = true,
    remotes = nothing,
)

if get(ENV, "CI", "") == "true"
    deploydocs(
        repo = "github.com/PenguinxCutCell/FrontTrackingMethods.jl",
        push_preview = true,
    )
end
