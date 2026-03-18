using FrontTrackingMethods
using StaticArrays

outdir = joinpath(@__DIR__, "output")
isdir(outdir) || mkpath(outdir)
outcsv = joinpath(outdir, "topology_merge_3d_study.csv")

patch_h_vals = [0.8, 1.0, 1.2]
d_trigger_vals = [1.8, 2.2, 2.6]
scopes = [:whole_component]
remesh_flags = [false]

open(outcsv, "w") do io
    println(io, "patch_h_factor,d_trigger,reconstruct_scope,remesh_after,changed,event,n_before,n_after,volume_drift,runtime_s")

    for hf in patch_h_vals, dt in d_trigger_vals, scope in scopes, remesh_after in remesh_flags
        fronts = make_two_spheres_merge_setup(c1=SVector(0.44, 0.5, 0.5), c2=SVector(0.56, 0.5, 0.5), R=0.10, refinement=1)
        state = MultiFrontState(fronts)
        V0 = front_enclosed_measure(state)

        handler = LocalCartesianTopologyHandler(
            d_trigger=dt,
            d_merge=dt,
            d_split=0.2,
            patch_h_factor=hf,
            reconstruct_scope=scope,
            max_components_created=8,
        )

        t0 = time()
        report = handle_topology_change!(state, handler)
        rt = time() - t0

        V1 = front_enclosed_measure(state)
        drift = abs(V1 - V0) / max(abs(V0), eps(Float64))
        println(io, "$hf,$dt,$scope,$remesh_after,$(report.changed),$(report.event_type),2,$(ncomponents(state)),$drift,$rt")
    end
end

println("Wrote: $outcsv")
