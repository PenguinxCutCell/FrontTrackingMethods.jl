# plotting_stubs.jl – Optional plotting API fallbacks.

const _MAKIE_ERR = "Makie plotting requires loading a Makie backend, e.g. `using CairoMakie`."

makie_theme(args...; kwargs...) = error(_MAKIE_ERR)
set_makie_theme!(args...; kwargs...) = error(_MAKIE_ERR)
plot_front(args...; kwargs...) = error(_MAKIE_ERR)
plot_state(args...; kwargs...) = error(_MAKIE_ERR)
plot_equation(args...; kwargs...) = error(_MAKIE_ERR)
animate_equation!(args...; kwargs...) = error(_MAKIE_ERR)
record_evolution!(args...; kwargs...) = error(_MAKIE_ERR)
snapshot(args...; kwargs...) = error(_MAKIE_ERR)
