using FrontTrackingMethods
using StaticArrays

outdir = joinpath(@__DIR__, "output")
isdir(outdir) || mkpath(outdir)
outcsv = joinpath(outdir, "topology_merge_2d_study.csv")

patch_h_vals = [0.4, 0.6, 0.8]
d_trigger_vals = [1.5, 2.0, 2.5]
scopes = [:local_patch, :whole_component]
remesh_flags = [false, true]

open(outcsv, "w") do io
    println(io, "patch_h_factor,d_trigger,reconstruct_scope,remesh_after,changed,event,n_before,n_after,area_drift,runtime_s")

    for hf in patch_h_vals, dt in d_trigger_vals, scope in scopes, remesh_after in remesh_flags
        fronts = make_two_circles_merge_setup(c1=SVector(0.47, 0.5), c2=SVector(0.53, 0.5), R=0.08, N=128)
        state = MultiFrontState(fronts)
        A0 = front_enclosed_measure(state)

        handler = LocalCartesianTopologyHandler(
            d_trigger=dt,
            d_merge=dt,
            d_split=0.8,
            patch_h_factor=hf,
            reconstruct_scope=scope,
            max_components_created=8,
        )

        t0 = time()
        report = handle_topology_change!(state, handler)
        if remesh_after
            for comp in eachcomponent(state)
                cstate = FrontState(comp.mesh; t=state.t)
                redistribute!(cstate, CurveEqualArcRedistributor())
                comp.mesh = cstate.mesh
                comp.geom = cstate.geom
            end
            refresh_geometry!(state; which=:all)
        end
        rt = time() - t0

        A1 = front_enclosed_measure(state)
        drift = abs(A1 - A0) / max(abs(A0), eps(Float64))
        println(io, "$hf,$dt,$scope,$remesh_after,$(report.changed),$(report.event_type),2,$(ncomponents(state)),$drift,$rt")
    end
end

println("Wrote: $outcsv")
