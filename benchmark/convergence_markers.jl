"""
Marker-count convergence study for all benchmark cases.

Runs each case over a range of N values, collects errors (E_area, E_shape, E_sym)
vs. mean edge length h, estimates log-log convergence orders via least-squares,
prints a formatted table, saves a CSV, and (if save_visuals=true) saves
log-log convergence plots to benchmark/output/convergence/.

Usage (from repo root):
    julia --project=. benchmark/convergence_markers.jl

Optional keyword arguments to run_marker_convergence():
    resolutions   – vector of N values (default: [32, 64, 128, 256, 512])
    save_visuals  – whether to save convergence plots (default: true)
"""

using DelimitedFiles
using Printf
using Statistics

include("translation_circle.jl")
include("solid_body_rotation_circle.jl")
include("zalesak_disk_rotation.jl")
include("single_vortex_circle.jl")

# ---------------------------------------------------------------------------
# Convergence rate estimator (least-squares log-log slope)
# ---------------------------------------------------------------------------

"""
    convergence_order(hs, errs) -> Float64

Estimate convergence order p such that E ≈ C * h^p via a least-squares fit
in log-log space.  Returns NaN if fewer than 2 valid (positive) data points.
"""
function convergence_order(hs::AbstractVector{<:Real}, errs::AbstractVector{<:Real})
    valid = [(h, e) for (h, e) in zip(hs, errs) if h > 0 && e > 0]
    length(valid) < 2 && return NaN
    logH = log.(first.(valid))
    logE = log.(last.(valid))
    # slope = cov(logH, logE) / var(logH)
    mH = mean(logH); mE = mean(logE)
    num = sum((lh - mH) * (le - mE) for (lh, le) in zip(logH, logE))
    den = sum((lh - mH)^2 for lh in logH)
    den ≈ 0 && return NaN
    return num / den
end

# ---------------------------------------------------------------------------
# Formatted table printer
# ---------------------------------------------------------------------------

function _print_convergence_table(case, rows; case_label="N")
    col_w = 10
    println()
    println("=" ^ 78)
    println("  Case: $case")
    println("=" ^ 78)
    hdr = @sprintf("  %-6s  %-12s  %-14s  %-14s  %-14s", case_label, "h", "E_area", "E_shape", "E_sym")
    println(hdr)
    println("-" ^ 78)
    for r in rows
        N_val  = r.N_val
        h      = r.h
        E_area = r.E_area
        E_shape = r.E_shape
        E_sym  = r.E_sym
        s = @sprintf("  %-6d  %-12.5e  %-14.6e  %-14.6e  %-14.6e",
                     N_val, h, E_area,
                     isnan(E_shape) ? NaN : E_shape,
                     E_sym)
        println(s)
    end
    println("-" ^ 78)
    hs     = [r.h      for r in rows]
    areas  = [r.E_area  for r in rows]
    shapes = [r.E_shape for r in rows]
    syms   = [r.E_sym   for r in rows]
    p_area  = convergence_order(hs, areas)
    p_shape = convergence_order(hs, shapes)
    p_sym   = convergence_order(hs, syms)
    @printf("  Convergence order:              p_area=%-6.2f  p_shape=%-6.2f  p_sym=%-6.2f\n",
            p_area, p_shape, p_sym)
    println("=" ^ 78)
end

# ---------------------------------------------------------------------------
# Per-case convergence sweeps
# ---------------------------------------------------------------------------

function convergence_translation(; resolutions, save_visuals::Bool)
    rows = []
    for N in resolutions
        r = run_translation_circle(N=N, save_visuals=false)
        push!(rows, (N_val=N, h=r.h, E_area=r.E_area, E_shape=r.E_shape, E_sym=r.E_sym))
    end
    _print_convergence_table("translation_circle (reversal)", rows)
    return rows
end

function convergence_solid_body(; resolutions, save_visuals::Bool)
    rows = []
    for N in resolutions
        r = run_solid_body_rotation_circle(N=N, save_visuals=false)
        push!(rows, (N_val=N, h=r.h, E_area=r.E_area, E_shape=r.E_shape, E_sym=r.E_sym))
    end
    _print_convergence_table("solid_body_rotation_circle", rows)
    return rows
end

function convergence_zalesak(; resolutions, save_visuals::Bool)
    rows = []
    for N in resolutions
        r = run_zalesak_disk_rotation(N_arc=N, save_visuals=false)
        push!(rows, (N_val=N, h=r.h, E_area=r.E_area, E_shape=NaN, E_sym=r.E_sym))
    end
    _print_convergence_table("zalesak_disk_rotation", rows; case_label="N_arc")
    return rows
end

function convergence_single_vortex(; resolutions, save_visuals::Bool)
    rows = []
    for N in resolutions
        r = run_single_vortex_circle(N=N, save_visuals=false)
        push!(rows, (N_val=N, h=r.h, E_area=r.E_area, E_shape=r.E_shape, E_sym=r.E_sym))
    end
    _print_convergence_table("single_vortex_circle (Rider-Kothe)", rows)
    return rows
end

# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------

function _write_convergence_csv(outfile, all_rows)
    open(outfile, "w") do io
        println(io, "case,N,h,E_area,E_shape,E_sym,p_area,p_shape,p_sym")
        for (case_name, rows) in all_rows
            hs     = [r.h       for r in rows]
            areas  = [r.E_area  for r in rows]
            shapes = [r.E_shape for r in rows]
            syms   = [r.E_sym   for r in rows]
            p_area  = convergence_order(hs, areas)
            p_shape = convergence_order(hs, shapes)
            p_sym   = convergence_order(hs, syms)
            for r in rows
                @printf(io, "%s,%d,%.8e,%.8e,%.8e,%.8e,%.4f,%.4f,%.4f\n",
                        case_name, r.N_val, r.h,
                        r.E_area, r.E_shape, r.E_sym,
                        p_area, p_shape, p_sym)
            end
        end
    end
    println("Saved convergence CSV: $outfile")
end

# ---------------------------------------------------------------------------
# Convergence plots (optional, requires CairoMakie in the environment)
# ---------------------------------------------------------------------------

function _plot_convergence_case(case_name, rows, outdir)
    # CairoMakie is already loaded via visualization.jl (included through benchmark drivers above)
    CairoMakie.activate!()

    hs     = [r.h       for r in rows]
    areas  = [r.E_area  for r in rows]
    shapes = [r.E_shape for r in rows]
    syms   = [r.E_sym   for r in rows]

    # Reference slopes: h^1, h^2
    hrange = [minimum(hs), maximum(hs)]
    ref1  = hrange ./ hrange[end] .* maximum(filter(isfinite, [areas; shapes; syms]))
    ref2  = (hrange ./ hrange[end]).^2 .* maximum(filter(isfinite, [areas; shapes; syms]))

    fig = Figure(size=(680, 500))
    ax  = Axis(fig[1, 1];
        title  = "Marker convergence: $case_name",
        xlabel = "mean edge length h",
        ylabel = "error",
        xscale = log10,
        yscale = log10,
        xticks = LogTicks(LinearTicks(6)),
        yticks = LogTicks(LinearTicks(6)),
    )

    # Reference lines
    lines!(ax, hrange, ref1; linestyle=:dash, color=:gray40, label="∝ h¹")
    lines!(ax, hrange, ref2; linestyle=:dot,  color=:gray70, label="∝ h²")

    # Error curves
    any(e -> e > 0, areas)  && scatterlines!(ax, hs, areas;  color=:royalblue,   marker=:circle,   label="E_area")
    any(e -> e > 0, syms)   && scatterlines!(ax, hs, syms;   color=:orangered,   marker=:utriangle,label="E_sym")
    valid_shapes = filter(isfinite, shapes)
    if !isempty(valid_shapes) && any(e -> e > 0, valid_shapes)
        scatterlines!(ax, hs, shapes; color=:green4, marker=:rect, label="E_shape")
    end

    axislegend(ax; position=:rb)
    save(joinpath(outdir, "$(case_name)_convergence.png"), fig)
    return fig
end

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

function run_marker_convergence(;
    resolutions::AbstractVector{Int} = [32, 64, 128, 256, 512],
    save_visuals::Bool = true,
)
    println()
    println("╔══════════════════════════════════════════════════════════════════════════╗")
    println("║         Marker-count convergence study — FrontTrackingMethods           ║")
    println("╚══════════════════════════════════════════════════════════════════════════╝")
    @printf("  Resolutions: %s\n", join(string.(resolutions), ", "))

    rows_trans  = convergence_translation(;  resolutions, save_visuals)
    rows_solid  = convergence_solid_body(;   resolutions, save_visuals)
    rows_zales  = convergence_zalesak(;      resolutions, save_visuals)
    rows_vortex = convergence_single_vortex(; resolutions, save_visuals)

    outdir = joinpath(@__DIR__, "output", "convergence")
    isdir(outdir) || mkpath(outdir)

    all_rows = [
        ("translation_circle",        rows_trans),
        ("solid_body_rotation_circle", rows_solid),
        ("zalesak_disk_rotation",      rows_zales),
        ("single_vortex_circle",       rows_vortex),
    ]

    outfile = joinpath(outdir, "marker_convergence.csv")
    _write_convergence_csv(outfile, all_rows)

    if save_visuals
        try
            for (name, rows) in all_rows
                _plot_convergence_case(name, rows, outdir)
            end
            println("Saved convergence plots to: $outdir")
        catch e
            @warn "Could not save convergence plots (CairoMakie may not be available)" exception=e
        end
    end

    return all_rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_marker_convergence()
end
