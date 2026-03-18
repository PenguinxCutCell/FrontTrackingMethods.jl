module MakieExt

using Makie
import FrontTrackingMethods as FTM
import FrontIntrinsicOps as FIO

function __init__()
    @info "FrontTrackingMethods Makie extension loaded"
end

_fio_has_makie_ext() = Base.get_extension(FIO, :MakieExt) !== nothing

# -----------------------------------------------------------------------------
# Theme
# -----------------------------------------------------------------------------

function FTM.makie_theme()
    if _fio_has_makie_ext()
        return FIO.makie_theme()
    end
    return Makie.Theme(
        fontsize = 14,
        backgroundcolor = :white,
        Axis = (
            xlabel = "x",
            ylabel = "y",
            aspect = Makie.DataAspect(),
            xgridvisible = true,
            ygridvisible = true,
        ),
        Axis3 = (
            xlabel = "x",
            ylabel = "y",
            zlabel = "z",
        ),
    )
end

function FTM.set_makie_theme!()
    Makie.set_theme!(FTM.makie_theme())
    return nothing
end

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

_curve_points_for_lines(mesh::FIO.CurveMesh) = begin
    pts = [Makie.Point2f(p[1], p[2]) for p in mesh.points]
    isempty(pts) || push!(pts, first(pts))
    pts
end

_surface_geometrybasics_mesh(mesh::FIO.SurfaceMesh) = begin
    pts = [FIO.GeometryBasics.Point3f(Float32(p[1]), Float32(p[2]), Float32(p[3])) for p in mesh.points]
    faces = [FIO.GeometryBasics.TriangleFace{Int}(f[1], f[2], f[3]) for f in mesh.faces]
    FIO.GeometryBasics.Mesh(pts, faces)
end

function _clear_axis_keep_decorations!(ax)
    for plt in copy(ax.scene.plots)
        delete!(ax.scene, plt)
    end
    return ax
end

function _apply_axis_limits!(ax; xlims=nothing, ylims=nothing, zlims=nothing)
    xlims === nothing || Makie.xlims!(ax, xlims...)
    ylims === nothing || Makie.ylims!(ax, ylims...)
    if zlims !== nothing && ax isa Makie.Axis3
        Makie.zlims!(ax, zlims...)
    end
    return ax
end

function _figure_axis_2d(; figure=nothing, axis=nothing, title=nothing)
    if axis !== nothing
        if title !== nothing
            axis.title = title
        end
        return figure === nothing ? axis.figure : figure, axis
    end
    fig = figure === nothing ? Makie.Figure() : figure
    ax = Makie.Axis(fig[1, 1]; title=title === nothing ? "" : title)
    return fig, ax
end

function _figure_axis_3d(; figure=nothing, axis=nothing, title=nothing)
    if axis !== nothing
        if title !== nothing
            axis.title = title
        end
        return figure === nothing ? axis.figure : figure, axis
    end
    fig = figure === nothing ? Makie.Figure() : figure
    ax = Makie.Axis3(fig[1, 1]; title=title === nothing ? "" : title, xlabel="x", ylabel="y", zlabel="z")
    return fig, ax
end

function _plot_curve_direct(mesh::FIO.CurveMesh;
    figure=nothing,
    axis=nothing,
    clear_axis::Bool=false,
    title=nothing,
    color=:royalblue,
    show_vertices::Bool=false,
    kwargs...,
)
    fig, ax = _figure_axis_2d(; figure=figure, axis=axis, title=title)
    clear_axis && _clear_axis_keep_decorations!(ax)
    pts = _curve_points_for_lines(mesh)
    p = Makie.lines!(ax, pts; color=color, kwargs...)
    if show_vertices
        Makie.scatter!(ax, [Makie.Point2f(q[1], q[2]) for q in mesh.points]; color=:black, markersize=6)
    end
    return fig, ax, p
end

function _plot_surface_direct(mesh::FIO.SurfaceMesh;
    figure=nothing,
    axis=nothing,
    clear_axis::Bool=false,
    title=nothing,
    color=:royalblue,
    wireframe::Bool=false,
    show_vertices::Bool=false,
    kwargs...,
)
    fig, ax = _figure_axis_3d(; figure=figure, axis=axis, title=title)
    clear_axis && _clear_axis_keep_decorations!(ax)
    gb = _surface_geometrybasics_mesh(mesh)
    p = Makie.mesh!(ax, gb; color=color, kwargs...)
    if wireframe
        Makie.wireframe!(ax, gb; color=:black, linewidth=1)
    end
    if show_vertices
        pts = [Makie.Point3f(q[1], q[2], q[3]) for q in mesh.points]
        Makie.scatter!(ax, pts; color=:black, markersize=4)
    end
    return fig, ax, p
end

function _vertex_scalar_values(state::FTM.FrontState, field)
    vals = if field === nothing
        nothing
    elseif field isa Symbol
        f = FTM.get_field(state, field)
        FTM.location(f) === :vertex || error("Only :vertex FrontField plotting is currently supported.")
        Base.values(f)
    elseif field isa FTM.FrontField
        FTM.location(field) === :vertex || error("Only :vertex FrontField plotting is currently supported.")
        Base.values(field)
    elseif field isa AbstractVector
        field
    else
        error("Unsupported field input type $(typeof(field)). Use Symbol, FrontField, or AbstractVector.")
    end

    vals === nothing && return nothing
    length(vals) == length(state.mesh.points) || error("Field length does not match number of mesh vertices.")
    eltype(vals) <: Number || error("Only scalar vertex fields are currently supported for coloring.")
    return vals
end

# -----------------------------------------------------------------------------
# Public plotting helpers
# -----------------------------------------------------------------------------

FTM.plot_front(mesh::Union{FIO.CurveMesh,FIO.SurfaceMesh}; kwargs...) = FIO.plot_front(mesh; kwargs...)
FTM.plot_front(state::FTM.FrontState; kwargs...) = FTM.plot_state(state; kwargs...)
FTM.plot_front(state::FTM.MultiFrontState; kwargs...) = FTM.plot_state(state; kwargs...)
FTM.plot_front(eq::FTM.FrontEquation; kwargs...) = FTM.plot_equation(eq; kwargs...)

function FTM.plot_front(mesh::FIO.CurveMesh;
    figure=nothing,
    axis=nothing,
    clear_axis::Bool=false,
    title=nothing,
    color=:royalblue,
    show_vertices::Bool=false,
    xlims=nothing,
    ylims=nothing,
    kwargs...,
)
    fig, ax, p = if _fio_has_makie_ext()
        FIO.plot_front(mesh;
            figure=figure, axis=axis, clear_axis=clear_axis, title=title,
            color=color, show_vertices=show_vertices, kwargs...)
    else
        _plot_curve_direct(mesh;
            figure=figure, axis=axis, clear_axis=clear_axis, title=title,
            color=color, show_vertices=show_vertices, kwargs...)
    end
    _apply_axis_limits!(ax; xlims=xlims, ylims=ylims)
    return fig, ax, p
end

function FTM.plot_front(mesh::FIO.SurfaceMesh;
    figure=nothing,
    axis=nothing,
    clear_axis::Bool=false,
    title=nothing,
    color=:royalblue,
    wireframe::Bool=false,
    show_vertices::Bool=false,
    xlims=nothing,
    ylims=nothing,
    zlims=nothing,
    kwargs...,
)
    fig, ax, p = if _fio_has_makie_ext()
        FIO.plot_front(mesh;
            figure=figure, axis=axis, clear_axis=clear_axis, title=title,
            color=color, wireframe=wireframe, show_vertices=show_vertices, kwargs...)
    else
        _plot_surface_direct(mesh;
            figure=figure, axis=axis, clear_axis=clear_axis, title=title,
            color=color, wireframe=wireframe, show_vertices=show_vertices, kwargs...)
    end
    _apply_axis_limits!(ax; xlims=xlims, ylims=ylims, zlims=zlims)
    return fig, ax, p
end

function FTM.plot_state(state::FTM.FrontState;
    figure=nothing,
    axis=nothing,
    clear_axis::Bool=false,
    show_vertices::Bool=false,
    show_normals::Bool=false,
    normal_scale::Real=0.05,
    normal_every::Int=1,
    field=nothing,
    wireframe::Bool=false,
    title=nothing,
    color=:royalblue,
    xlims=nothing,
    ylims=nothing,
    zlims=nothing,
    kwargs...,
)
    vals = _vertex_scalar_values(state, field)
    mesh = state.mesh

    if mesh isa FIO.CurveMesh
        fig, ax, p = FTM.plot_front(mesh;
            figure=figure,
            axis=axis,
            clear_axis=clear_axis,
            title=title,
            color=color,
            show_vertices=show_vertices,
            xlims=xlims,
            ylims=ylims,
            kwargs...,
        )
        if vals !== nothing
            pts = [Makie.Point2f(q[1], q[2]) for q in mesh.points]
            Makie.scatter!(ax, pts; color=vals, colormap=:viridis, markersize=7)
        end
        if show_normals
            if _fio_has_makie_ext()
                FIO.plot_normals(mesh, state.geom; figure=fig, axis=ax, scale=normal_scale, every=normal_every)
            end
        end
        _apply_axis_limits!(ax; xlims=xlims, ylims=ylims)
        return fig, ax, p
    end

    # Surface case
    surface_color = vals === nothing ? color : vals
    fig, ax, p = FTM.plot_front(mesh;
        figure=figure,
        axis=axis,
        clear_axis=clear_axis,
        title=title,
        color=surface_color,
        wireframe=wireframe,
        show_vertices=show_vertices,
        xlims=xlims,
        ylims=ylims,
        zlims=zlims,
        kwargs...,
    )
    if show_normals
        if _fio_has_makie_ext()
            FIO.plot_normals(mesh, state.geom; figure=fig, axis=ax, scale=normal_scale, every=normal_every)
        end
    end
    _apply_axis_limits!(ax; xlims=xlims, ylims=ylims, zlims=zlims)
    return fig, ax, p
end

function FTM.plot_state(state::FTM.MultiFrontState;
    figure=nothing,
    axis=nothing,
    clear_axis::Bool=false,
    show_vertices::Bool=false,
    show_normals::Bool=false,
    normal_scale::Real=0.05,
    normal_every::Int=1,
    field=nothing,
    wireframe::Bool=false,
    title=nothing,
    color=:royalblue,
    xlims=nothing,
    ylims=nothing,
    zlims=nothing,
    kwargs...,
)
    FTM.ncomponents(state) > 0 || error("plot_state(::MultiFrontState): state has no components.")

    comps = collect(FTM.eachcomponent(state))
    first_mesh = comps[1].mesh
    is2d = first_mesh isa FIO.CurveMesh

    fig, ax = if is2d
        _figure_axis_2d(; figure=figure, axis=axis, title=title)
    else
        _figure_axis_3d(; figure=figure, axis=axis, title=title)
    end
    clear_axis && _clear_axis_keep_decorations!(ax)

    last_plot = nothing
    for comp in comps
        cstate = FTM.FrontState{typeof(comp.mesh),typeof(comp.geom),typeof(comp.dec)}(
            comp.mesh,
            comp.geom,
            comp.dec,
            FTM.current_time(state),
            comp.fields,
            comp.cache,
        )

        _, _, p = FTM.plot_state(cstate;
            figure=fig,
            axis=ax,
            clear_axis=false,
            show_vertices=show_vertices,
            show_normals=show_normals,
            normal_scale=normal_scale,
            normal_every=normal_every,
            field=field,
            wireframe=wireframe,
            title=nothing,
            color=color,
            xlims=xlims,
            ylims=ylims,
            zlims=zlims,
            kwargs...,
        )
        last_plot = p
    end

    if title !== nothing
        ax.title = title
    end
    _apply_axis_limits!(ax; xlims=xlims, ylims=ylims, zlims=zlims)
    return fig, ax, last_plot
end

function FTM.plot_equation(eq::FTM.FrontEquation; title=nothing, kwargs...)
    tnow = FTM.current_time(eq)
    tt = title === nothing ? "t = $(round(tnow, digits=4))" : title
    return FTM.plot_state(FTM.current_state(eq); title=tt, kwargs...)
end

function FTM.snapshot(obj, filename::AbstractString; kwargs...)
    fig = if obj isa FIO.CurveMesh || obj isa FIO.SurfaceMesh
        f, _, _ = FTM.plot_front(obj; kwargs...)
        f
    elseif obj isa FTM.FrontState
        f, _, _ = FTM.plot_state(obj; kwargs...)
        f
    elseif obj isa FTM.MultiFrontState
        f, _, _ = FTM.plot_state(obj; kwargs...)
        f
    elseif obj isa FTM.FrontEquation
        f, _, _ = FTM.plot_equation(obj; kwargs...)
        f
    else
        error("snapshot: unsupported object type $(typeof(obj)).")
    end
    Makie.save(filename, fig)
    return filename
end

function FTM.record_evolution!(eq::FTM.FrontEquation, filename::AbstractString, times;
    show_vertices::Bool=false,
    show_normals::Bool=false,
    normal_scale::Real=0.05,
    normal_every::Int=1,
    field=nothing,
    wireframe::Bool=false,
    color=:royalblue,
    xlims=nothing,
    ylims=nothing,
    zlims=nothing,
    kwargs...,
)
    isempty(times) && error("record_evolution!: `times` must be non-empty.")

    fig, ax, _ = FTM.plot_equation(eq;
        show_vertices=show_vertices,
        show_normals=show_normals,
        normal_scale=normal_scale,
        normal_every=normal_every,
        field=field,
        wireframe=wireframe,
        color=color,
        xlims=xlims,
        ylims=ylims,
        zlims=zlims,
        kwargs...,
    )

    Makie.record(fig, filename, collect(times)) do ttarget
        FTM.integrate!(eq, Float64(ttarget))
        _clear_axis_keep_decorations!(ax)
        FTM.plot_state(FTM.current_state(eq);
            figure=fig,
            axis=ax,
            clear_axis=false,
            title="t = $(round(FTM.current_time(eq), digits=4))",
            show_vertices=show_vertices,
            show_normals=show_normals,
            normal_scale=normal_scale,
            normal_every=normal_every,
            field=field,
            wireframe=wireframe,
            color=color,
            xlims=xlims,
            ylims=ylims,
            zlims=zlims,
            kwargs...,
        )
    end

    return filename
end

FTM.animate_equation!(eq::FTM.FrontEquation, filename::AbstractString, times; kwargs...) =
    FTM.record_evolution!(eq, filename, times; kwargs...)

# Convenience user API
Makie.plot(mesh::Union{FIO.CurveMesh,FIO.SurfaceMesh}; kwargs...) = first(FTM.plot_front(mesh; kwargs...))
Makie.plot(state::FTM.FrontState; kwargs...) = first(FTM.plot_state(state; kwargs...))
Makie.plot(state::FTM.MultiFrontState; kwargs...) = first(FTM.plot_state(state; kwargs...))
Makie.plot(eq::FTM.FrontEquation; kwargs...) = first(FTM.plot_equation(eq; kwargs...))

end
