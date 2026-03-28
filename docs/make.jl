using Documenter
using FrontTrackingMethods

function copy_tree!(src::AbstractString, dst::AbstractString)
    mkpath(dst)
    for (root, dirs, files) in walkdir(src)
        rel = relpath(root, src)
        out_root = rel == "." ? dst : joinpath(dst, rel)
        mkpath(out_root)
        for d in dirs
            mkpath(joinpath(out_root, d))
        end
        for f in files
            cp(joinpath(root, f), joinpath(out_root, f); force=true)
        end
    end
    return nothing
end

function sync_output_assets!()
    repo_root = normpath(joinpath(@__DIR__, ".."))
    docs_src = joinpath(@__DIR__, "src")
    generated_root = joinpath(docs_src, "assets", "generated")

    mkpath(generated_root)

    output_pairs = [
        (joinpath(repo_root, "examples", "output"), joinpath(generated_root, "examples_output")),
        (joinpath(repo_root, "benchmark", "output"), joinpath(generated_root, "benchmark_output")),
    ]

    for (src, dst) in output_pairs
        isdir(dst) && rm(dst; recursive=true, force=true)
        if isdir(src)
            copy_tree!(src, dst)
        else
            @warn "Output folder not found; skipping docs asset sync" src
        end
    end

    return nothing
end

sync_output_assets!()

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
