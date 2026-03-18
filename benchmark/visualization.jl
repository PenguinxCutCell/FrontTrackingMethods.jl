using CairoMakie

function benchmark_output_dir(case_name::String)
    outdir = joinpath(@__DIR__, "output", case_name)
    isdir(outdir) || mkpath(outdir)
    return outdir
end

function default_times(Tfinal::Float64; nframes::Int=160)
    return collect(range(0.0, stop=Tfinal, length=nframes))
end

function save_case_visuals!(
    case_name::String,
    eq_anim,
    Tfinal::Float64;
    nframes::Int=160,
    xlims=nothing,
    ylims=nothing,
    zlims=nothing,
)
    CairoMakie.activate!()

    FrontTrackingMethods.set_makie_theme!()
    outdir = benchmark_output_dir(case_name)

    snapshot(eq_anim, joinpath(outdir, "initial.png");
        title="$case_name: initial",
        xlims=xlims,
        ylims=ylims,
        zlims=zlims,
    )
    record_evolution!(eq_anim, joinpath(outdir, "animation.mp4"), default_times(Tfinal; nframes=nframes);
        title=case_name,
        xlims=xlims,
        ylims=ylims,
        zlims=zlims,
    )
    snapshot(eq_anim, joinpath(outdir, "final.png");
        title="$case_name: final",
        xlims=xlims,
        ylims=ylims,
        zlims=zlims,
    )
    println("Saved visuals to: $outdir")
    return outdir
end
