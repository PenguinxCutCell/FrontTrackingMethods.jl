using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using Printf

function _surface_area(mesh)
    FrontIntrinsicOps.measure(mesh, FrontIntrinsicOps.compute_geometry(mesh))
end

function _run_enright_case(mesh0, u, tf, dt; label, redistributor=nothing, log_every=2)
    eq = FrontEquation(;
        terms=AdvectionTerm(u),
        front=deepcopy(mesh0),
        integrator=RK2(),
        redistributor=redistributor,
    )

    min_angle_hist = Float64[]
    degen_hist = Float64[]
    t_hist = Float64[]

    cb = EveryNSteps(log_every, (state, t, step) -> begin
        q = surface_quality_summary(state.mesh; degenerate_atol=1e-12)
        push!(min_angle_hist, q.angle.min_angle)
        push!(degen_hist, q.degenerate_fraction)
        push!(t_hist, t)
        nothing
    end)

    t0 = time()
    integrate!(eq, tf; dt=dt, callback=cb)
    runtime = time() - t0

    state = current_state(eq)
    meshf = state.mesh

    V0 = enclosed_measure(mesh0)
    Vf = enclosed_measure(meshf)
    A0 = _surface_area(mesh0)
    Af = _surface_area(meshf)
    dH = symmetric_hausdorff_surface(meshf, mesh0)
    qf = surface_quality_summary(meshf; degenerate_atol=1e-12)

    return (
        label=label,
        runtime_s=runtime,
        volume_rel=abs(Vf - V0) / max(abs(V0), eps()),
        area_rel=abs(Af - A0) / max(abs(A0), eps()),
        front_error=dH,
        min_angle_final=qf.angle.min_angle,
        min_angle_worst=isempty(min_angle_hist) ? qf.angle.min_angle : minimum(min_angle_hist),
        degenerate_final=qf.degenerate_fraction,
        degenerate_worst=isempty(degen_hist) ? qf.degenerate_fraction : maximum(degen_hist),
        n_samples=length(t_hist),
        split_count=NaN,
        collapse_count=NaN,
        flip_count=NaN,
    )
end

function run_enright3d_remeshing_study(; refinement=2, tf=1.5, dt=0.05, write_csv::Bool=true)
    mesh0 = make_sphere_benchmark_surface(center=SVector(0.35, 0.35, 0.35), R=0.15, refinement=refinement)
    u = enright_3d(; T=3.0)

    cases = [
        ("none", nothing),
        ("tangential_i2", SurfaceTangentialRedistributor(; iterations=2, strength=0.35)),
        ("experimental_i1", ExperimentalSurfaceRemesher(; iterations=1, strength=0.2, volume_correction=true)),
        ("experimental_i2", ExperimentalSurfaceRemesher(; iterations=2, strength=0.25, volume_correction=true)),
    ]

    rows = NamedTuple[]
    for (label, redist) in cases
        push!(rows, _run_enright_case(mesh0, u, tf, dt; label=label, redistributor=redist, log_every=2))
    end

    println("case,runtime_s,volume_rel,area_rel,front_error,min_angle_final,min_angle_worst,degenerate_final,degenerate_worst,n_samples,split_count,collapse_count,flip_count")
    for r in rows
        @printf("%s,%.4f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%d,NaN,NaN,NaN\n",
            r.label, r.runtime_s, r.volume_rel, r.area_rel, r.front_error,
            r.min_angle_final, r.min_angle_worst, r.degenerate_final, r.degenerate_worst, r.n_samples)
    end

    if write_csv
        outdir = joinpath(@__DIR__, "output")
        isdir(outdir) || mkpath(outdir)
        outfile = joinpath(outdir, "enright3d_remeshing_study.csv")
        open(outfile, "w") do io
            println(io, "case,runtime_s,volume_rel,area_rel,front_error,min_angle_final,min_angle_worst,degenerate_final,degenerate_worst,n_samples,split_count,collapse_count,flip_count")
            for r in rows
                @printf(io, "%s,%.8f,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%d,NaN,NaN,NaN\n",
                    r.label, r.runtime_s, r.volume_rel, r.area_rel, r.front_error,
                    r.min_angle_final, r.min_angle_worst, r.degenerate_final, r.degenerate_worst, r.n_samples)
            end
        end
        println("Saved: $outfile")
    end

    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_enright3d_remeshing_study()
end
