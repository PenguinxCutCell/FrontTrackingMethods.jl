using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Printf

function _surface_area(mesh)
    FrontIntrinsicOps.measure(mesh, FrontIntrinsicOps.compute_geometry(mesh))
end

function _slot_proxy(mesh_ref, mesh; center, x_tol=0.03, z_cut=0.08)
    pts = [p for p in mesh_ref.points if abs(p[1] - center[1]) <= x_tol && p[3] <= center[3] - z_cut]
    isempty(pts) && return NaN
    return maximum(minimum(norm(p - q) for q in mesh.points) for p in pts)
end

function run_zalesak_sphere_study(; refinements=(2, 3), dt=0.05, write_csv::Bool=true)
    rows = NamedTuple[]

    for refinement in refinements
        center_shape = SVector(0.5, 0.75, 0.5)
        mesh0 = make_zalesak_sphere_surface(
            center=center_shape,
            R=0.15,
            slot_width=0.05,
            slot_depth=0.125,
            refinement=refinement,
        )

        u = rigid_rotation_3d(; center=SVector(0.5, 0.5, 0.5), axis=:z, omega=2π)

        for (label, redist) in (
            ("none", nothing),
            ("experimental", ExperimentalSurfaceRemesher(; iterations=1, strength=0.15, volume_correction=false)),
        )
            eq = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2(), redistributor=redist)
            t0 = time()
            integrate!(eq, 1.0; dt=dt)
            runtime = time() - t0

            meshf = current_state(eq).mesh
            A0 = _surface_area(mesh0)
            Af = _surface_area(meshf)
            c0 = sum(mesh0.points) / length(mesh0.points)
            cf = sum(meshf.points) / length(meshf.points)
            dH = symmetric_hausdorff_surface(meshf, mesh0)
            q = surface_quality_summary(meshf; degenerate_atol=1e-12)
            slot_err = _slot_proxy(mesh0, meshf; center=center_shape)

            push!(rows, (
                refinement=refinement,
                remesher=label,
                runtime_s=runtime,
                closed=is_closed(mesh0),
                volume_rel=NaN,
                area_rel=abs(Af - A0) / max(abs(A0), eps()),
                centroid_err=norm(cf - c0),
                front_error=dH,
                slot_proxy=slot_err,
                min_angle=q.angle.min_angle,
                degenerate_fraction=q.degenerate_fraction,
            ))
        end
    end

    println("refinement,remesher,runtime_s,closed,volume_rel,area_rel,centroid_err,front_error,slot_proxy,min_angle,degenerate_fraction")
    for r in rows
        @printf("%d,%s,%.4f,%s,NaN,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
            r.refinement, r.remesher, r.runtime_s, string(r.closed), r.area_rel,
            r.centroid_err, r.front_error, r.slot_proxy, r.min_angle, r.degenerate_fraction)
    end

    if write_csv
        outdir = joinpath(@__DIR__, "output")
        isdir(outdir) || mkpath(outdir)
        outfile = joinpath(outdir, "zalesak_sphere_study.csv")
        open(outfile, "w") do io
            println(io, "refinement,remesher,runtime_s,closed,volume_rel,area_rel,centroid_err,front_error,slot_proxy,min_angle,degenerate_fraction")
            for r in rows
                @printf(io, "%d,%s,%.8f,%s,NaN,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                    r.refinement, r.remesher, r.runtime_s, string(r.closed), r.area_rel,
                    r.centroid_err, r.front_error, r.slot_proxy, r.min_angle, r.degenerate_fraction)
            end
        end
        println("Saved: $outfile")
    end

    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_zalesak_sphere_study()
end
