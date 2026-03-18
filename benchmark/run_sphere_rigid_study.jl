using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Printf

function _surface_area(mesh)
    FrontIntrinsicOps.measure(mesh, FrontIntrinsicOps.compute_geometry(mesh))
end

function _run_case(mesh0, u, tf, dt; redistributor=nothing)
    eq = FrontEquation(;
        terms=AdvectionTerm(u),
        front=deepcopy(mesh0),
        integrator=RK2(),
        redistributor=redistributor,
    )
    t0 = time()
    integrate!(eq, tf; dt=dt)
    runtime = time() - t0

    meshf = current_state(eq).mesh
    q = surface_quality_summary(meshf; degenerate_atol=1e-12)

    return (
        meshf=meshf,
        runtime=runtime,
        q=q,
    )
end

function run_sphere_rigid_study(; refinements=(2, 3), dt=0.05, write_csv::Bool=true)
    rows = NamedTuple[]

    for refinement in refinements
        mesh0 = make_sphere_benchmark_surface(center=SVector(0.35, 0.35, 0.35), R=0.15, refinement=refinement)

        # Translation case
        u_tr = rigid_translation_velocity(SVector(0.1, -0.05, 0.02))
        tf_tr = 1.0
        base_tr = _run_case(mesh0, u_tr, tf_tr, dt)
        rem_tr  = _run_case(mesh0, u_tr, tf_tr, dt;
            redistributor=ExperimentalSurfaceRemesher(; iterations=2, strength=0.2, volume_correction=true))

        for (label, out) in (("none", base_tr), ("experimental", rem_tr))
            meshf = out.meshf
            V0 = enclosed_measure(mesh0)
            Vf = enclosed_measure(meshf)
            A0 = _surface_area(mesh0)
            Af = _surface_area(meshf)
            dH = symmetric_hausdorff_surface(meshf, mesh0)
            c0 = sum(mesh0.points) / length(mesh0.points)
            cf = sum(meshf.points) / length(meshf.points)

            push!(rows, (
                study="translation",
                refinement=refinement,
                remesher=label,
                runtime_s=out.runtime,
                volume_rel=abs(Vf - V0) / max(abs(V0), eps()),
                area_rel=abs(Af - A0) / max(abs(A0), eps()),
                centroid_err=norm(cf - c0 - SVector(0.1, -0.05, 0.02) * tf_tr),
                front_error=dH,
                min_angle=out.q.angle.min_angle,
                degenerate_fraction=out.q.degenerate_fraction,
            ))
        end

        # Rotation case
        mesh_rot0 = make_sphere_benchmark_surface(center=SVector(0.5, 0.75, 0.5), R=0.15, refinement=refinement)
        u_rot = rigid_rotation_3d(; center=SVector(0.5, 0.5, 0.5), axis=:z, omega=2π)
        tf_rot = 1.0
        base_rot = _run_case(mesh_rot0, u_rot, tf_rot, dt)
        rem_rot  = _run_case(mesh_rot0, u_rot, tf_rot, dt;
            redistributor=ExperimentalSurfaceRemesher(; iterations=2, strength=0.2, volume_correction=true))

        for (label, out) in (("none", base_rot), ("experimental", rem_rot))
            meshf = out.meshf
            V0 = enclosed_measure(mesh_rot0)
            Vf = enclosed_measure(meshf)
            A0 = _surface_area(mesh_rot0)
            Af = _surface_area(meshf)
            dH = symmetric_hausdorff_surface(meshf, mesh_rot0)
            c0 = sum(mesh_rot0.points) / length(mesh_rot0.points)
            cf = sum(meshf.points) / length(meshf.points)

            push!(rows, (
                study="rotation",
                refinement=refinement,
                remesher=label,
                runtime_s=out.runtime,
                volume_rel=abs(Vf - V0) / max(abs(V0), eps()),
                area_rel=abs(Af - A0) / max(abs(A0), eps()),
                centroid_err=norm(cf - c0),
                front_error=dH,
                min_angle=out.q.angle.min_angle,
                degenerate_fraction=out.q.degenerate_fraction,
            ))
        end
    end

    println("study,refinement,remesher,runtime_s,volume_rel,area_rel,centroid_err,front_error,min_angle,degenerate_fraction")
    for r in rows
        @printf("%s,%d,%s,%.4f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
            r.study, r.refinement, r.remesher, r.runtime_s, r.volume_rel, r.area_rel,
            r.centroid_err, r.front_error, r.min_angle, r.degenerate_fraction)
    end

    if write_csv
        outdir = joinpath(@__DIR__, "output")
        isdir(outdir) || mkpath(outdir)
        outfile = joinpath(outdir, "sphere_rigid_study.csv")
        open(outfile, "w") do io
            println(io, "study,refinement,remesher,runtime_s,volume_rel,area_rel,centroid_err,front_error,min_angle,degenerate_fraction")
            for r in rows
                @printf(io, "%s,%d,%s,%.8f,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                    r.study, r.refinement, r.remesher, r.runtime_s, r.volume_rel, r.area_rel,
                    r.centroid_err, r.front_error, r.min_angle, r.degenerate_fraction)
            end
        end
        println("Saved: $outfile")
    end

    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_sphere_rigid_study()
end
