using FrontTrackingMethods
using StaticArrays
using Dates
using CairoMakie

example_name = "topology_merge_two_circles_2d"
outdir = joinpath(@__DIR__, "output", example_name)
isdir(outdir) || mkpath(outdir)

fronts = make_two_circles_merge_setup(
    c1=SVector(0.47, 0.5),
    c2=SVector(0.53, 0.5),
    R=0.08,
    N=128,
)
state = MultiFrontState(fronts)

A0 = front_enclosed_measure(state)

handler = LocalCartesianTopologyHandler(
    d_merge=2.2,
    d_split=0.8,
    patch_h_factor=0.45,
    patch_margin_factor=2.0,
    reconstruct_scope=:whole_component,
    preserve_fields=true,
)

report = handle_topology_change!(state, handler)
A1 = front_enclosed_measure(state)

println("[$(Dates.now())] $example_name")
println("components: 2 -> $(ncomponents(state))")
println("event: $(report.event_type), changed=$(report.changed)")
println("area drift: ", abs(A1 - A0) / max(abs(A0), eps(Float64)))

metrics_path = joinpath(outdir, "metrics.csv")
open(metrics_path, "w") do io
    println(io, "time,components,event,changed,area0,area1,area_drift")
    println(io, "0.0,2,none,false,$A0,$A0,0.0")
    println(io, "1.0,$(ncomponents(state)),$(report.event_type),$(report.changed),$A0,$A1,$(abs(A1-A0)/max(abs(A0),eps(Float64)))")
end

# Optional plotting if Makie backend is available.
        using CairoMakie
        FrontTrackingMethods.set_makie_theme!()

        eq0 = FrontEquation(terms=(AdvectionTerm(SVector(0.0, 0.0)),), state=MultiFrontState(fronts), topology_handler=NoTopologyChange())
        snapshot(eq0, joinpath(outdir, "initial.png"); title=example_name * " initial")

        eqf = FrontEquation(terms=(AdvectionTerm(SVector(0.0, 0.0)),), state=state, topology_handler=NoTopologyChange())
        snapshot(eqf, joinpath(outdir, "final.png"); title=example_name * " final")

