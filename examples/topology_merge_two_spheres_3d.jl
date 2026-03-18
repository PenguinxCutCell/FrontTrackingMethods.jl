using FrontTrackingMethods
using StaticArrays
using Dates

example_name = "topology_merge_two_spheres_3d"
outdir = joinpath(@__DIR__, "output", example_name)
isdir(outdir) || mkpath(outdir)

fronts = make_two_spheres_merge_setup(
    c1=SVector(0.44, 0.5, 0.5),
    c2=SVector(0.56, 0.5, 0.5),
    R=0.10,
    refinement=1,
)
state = MultiFrontState(fronts)

V0 = front_enclosed_measure(state)

handler = LocalCartesianTopologyHandler(
    d_merge=2.2,
    d_split=0.2,
    patch_h_factor=0.9,
    reconstruct_scope=:whole_component,
    max_components_created=8,
)

report = handle_topology_change!(state, handler)
V1 = front_enclosed_measure(state)

println("[$(Dates.now())] $example_name")
println("components: 2 -> $(ncomponents(state))")
println("event: $(report.event_type), changed=$(report.changed)")
println("volume drift: ", abs(V1 - V0) / max(abs(V0), eps(Float64)))

metrics_path = joinpath(outdir, "metrics.csv")
open(metrics_path, "w") do io
    println(io, "time,components,event,changed,volume0,volume1,volume_drift")
    println(io, "0.0,2,none,false,$V0,$V0,0.0")
    println(io, "1.0,$(ncomponents(state)),$(report.event_type),$(report.changed),$V0,$V1,$(abs(V1-V0)/max(abs(V0),eps(Float64)))")
end
