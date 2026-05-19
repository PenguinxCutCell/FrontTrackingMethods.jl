using DelimitedFiles
using LinearAlgebra
using Printf
using Statistics
using StaticArrays

using FrontIntrinsicOps
using FrontTrackingMethods

include("error_metrics.jl")

const BENCH_CENTER = SVector{2,Float64}(0.5, 0.75)
const ROT_CENTER = SVector{2,Float64}(0.5, 0.5)

function make_ellipse_benchmark_curve(;
    center::SVector{2,Float64}=BENCH_CENTER,
    a::Float64=0.30,
    b::Float64=0.15,
    N::Int=64,
)
    pts = [
        center + SVector{2,Float64}(a * cos(2π * k / N), b * sin(2π * k / N))
        for k in 0:N-1
    ]
    edges = [SVector{2,Int}(i, mod1(i + 1, N)) for i in 1:N]
    return CurveMesh{Float64}(pts, edges)
end

@inline function ellipse_curvature(a::Float64, b::Float64, θ::Float64)
    return (a * b) / (a^2 * sin(θ)^2 + b^2 * cos(θ)^2)^(1.5)
end

@inline function ellipse_inward_normal(p, center, a::Float64, b::Float64)
    x = p[1] - center[1]
    y = p[2] - center[2]
    return normalize(SVector{2,Float64}(-x / a^2, -y / b^2))
end

@inline function circle_inward_normal(p, center)
    return normalize(center - p)
end

function ramanujan_ellipse_perimeter(a::Float64, b::Float64)
    return π * (3 * (a + b) - sqrt((3a + b) * (a + 3b)))
end

function standard_curve_metrics(mesh::CurveMesh; exact_area::Float64, exact_perimeter::Float64)
    geom = compute_geometry(mesh)
    A_chord = enclosed_measure(mesh)
    A_arc = enclosed_measure(mesh; correction=:arc)
    P_chord = measure(mesh, geom)
    P_arc = measure(mesh, geom; correction=:arc)
    lengths = geom.edge_lengths
    return (
        E_A_shoe = abs(A_chord - exact_area) / exact_area,
        E_A_arc = abs(A_arc - exact_area) / exact_area,
        E_P_chord = abs(P_chord - exact_perimeter) / exact_perimeter,
        E_P_arc = abs(P_arc - exact_perimeter) / exact_perimeter,
        uniformity = std(lengths) / mean(lengths),
    )
end

function curvature_errors(mesh::CurveMesh, k_exact::AbstractVector{Float64})
    geom = compute_geometry(mesh)
    κ = FrontIntrinsicOps.menger_curvature(mesh)
    weights = geom.vertex_dual_lengths
    diff = κ .- k_exact
    denom2 = sqrt(sum(weights .* k_exact.^2))
    return (
        E_kappa_inf = maximum(abs.(diff) ./ max.(abs.(k_exact), eps(Float64))),
        E_kappa_2 = sqrt(sum(weights .* diff.^2)) / max(denom2, eps(Float64)),
    )
end

function max_normal_angle(mesh::CurveMesh, exact_normals::AbstractVector{SVector{2,Float64}}, method::Symbol)
    geom = compute_geometry(mesh)
    normals = if method === :laplacian
        dec = build_dec(mesh, geom)
        vertex_normals(mesh, geom; method=method, dec=dec)
    else
        vertex_normals(mesh, geom; method=method)
    end
    return maximum(acos(clamp(dot(normals[i], exact_normals[i]), -1.0, 1.0))
                   for i in eachindex(normals))
end

function ellipse_form_error(mesh::CurveMesh; center=BENCH_CENTER, a::Float64=0.30, b::Float64=0.15)
    return maximum(abs(((p[1] - center[1]) / a)^2 + ((p[2] - center[2]) / b)^2 - 1)
                   for p in mesh.points)
end

curve_vertex_centroid(mesh::CurveMesh) = sum(mesh.points) / length(mesh.points)

function run_static_circle_convergence(; Ns=(16, 32, 64, 128, 256), R::Float64=0.15)
    rows = NamedTuple[]
    exact_area = π * R^2
    exact_perimeter = 2π * R
    for N in Ns
        mesh = make_circle_benchmark_curve(center=BENCH_CENTER, R=R, N=N)
        exact_k = fill(1 / R, N)
        exact_normals = [circle_inward_normal(p, BENCH_CENTER) for p in mesh.points]
        m = standard_curve_metrics(mesh; exact_area=exact_area, exact_perimeter=exact_perimeter)
        ke = curvature_errors(mesh, exact_k)
        row = (
            case = "static_circle",
            N = N,
            m...,
            ke...,
            E_n_bisector = max_normal_angle(mesh, exact_normals, :bisector),
            E_n_osculating = max_normal_angle(mesh, exact_normals, :osculating),
            E_n_laplacian = max_normal_angle(mesh, exact_normals, :laplacian),
            E_form = maximum(abs(norm(p - BENCH_CENTER) - R) / R for p in mesh.points),
        )
        push!(rows, row)
    end
    return rows
end

function run_static_ellipse_convergence(; Ns=(32, 64, 128), a::Float64=0.30, b::Float64=0.15)
    rows = NamedTuple[]
    exact_area = π * a * b
    exact_perimeter = ramanujan_ellipse_perimeter(a, b)
    for N in Ns
        mesh = make_ellipse_benchmark_curve(a=a, b=b, N=N)
        θs = [2π * (i - 1) / N for i in 1:N]
        exact_k = [ellipse_curvature(a, b, θ) for θ in θs]
        exact_normals = [ellipse_inward_normal(p, BENCH_CENTER, a, b) for p in mesh.points]
        m = standard_curve_metrics(mesh; exact_area=exact_area, exact_perimeter=exact_perimeter)
        ke = curvature_errors(mesh, exact_k)
        row = (
            case = "static_ellipse",
            N = N,
            m...,
            ke...,
            E_n_bisector = max_normal_angle(mesh, exact_normals, :bisector),
            E_n_osculating = max_normal_angle(mesh, exact_normals, :osculating),
            E_n_laplacian = max_normal_angle(mesh, exact_normals, :laplacian),
            E_form = ellipse_form_error(mesh; a=a, b=b),
        )
        push!(rows, row)
    end
    return rows
end

function run_redistribution_circle_study(; Ns=(8, 16, 32, 64), R::Float64=0.15, iterations::Int=8)
    rows = NamedTuple[]
    exact_area = π * R^2
    exact_perimeter = 2π * R
    for N in Ns
        base = make_circle_benchmark_curve(center=BENCH_CENTER, R=R, N=N)
        pts = copy(base.points)
        for i in eachindex(pts)
            θ = 2π * (i - 1) / N
            pts[i] = BENCH_CENTER + (R * (1 + 0.12 * sin(3θ))) * SVector(cos(θ), sin(θ))
        end

        state = FrontState(CurveMesh{Float64}(pts, base.edges))
        before = standard_curve_metrics(state.mesh; exact_area=exact_area, exact_perimeter=exact_perimeter)
        redistribute!(state, PoissonTangentialRedistributor(; iterations=iterations, pseudo_dt=0.1, omega=1.0))
        after = standard_curve_metrics(state.mesh; exact_area=exact_area, exact_perimeter=exact_perimeter)

        push!(rows, (
            case = "poisson_redistribution_circle_before",
            N = N,
            before...,
            E_kappa_inf = NaN,
            E_kappa_2 = NaN,
            E_n_bisector = NaN,
            E_n_osculating = NaN,
            E_n_laplacian = NaN,
            E_form = maximum(abs(norm(p - BENCH_CENTER) - R) / R for p in pts),
        ))
        push!(rows, (
            case = "poisson_redistribution_circle_after",
            N = N,
            after...,
            E_kappa_inf = NaN,
            E_kappa_2 = NaN,
            E_n_bisector = NaN,
            E_n_osculating = NaN,
            E_n_laplacian = NaN,
            E_form = maximum(abs(norm(p - BENCH_CENTER) - R) / R for p in state.mesh.points),
        ))
    end
    return rows
end

function run_projected_rotation_circle_study(; Ns=(8, 16, 32, 64), R::Float64=0.15)
    rows = NamedTuple[]
    exact_area = π * R^2
    exact_perimeter = 2π * R
    u = rigid_rotation_2d(; center=ROT_CENTER, omega=2π)
    for N in Ns
        mesh0 = make_circle_benchmark_curve(center=BENCH_CENTER, R=R, N=N)
        h = mean_edge_size(mesh0)
        dt = 0.05 * h / π
        eq = FrontEquation(
            terms = ProjectedAdvectionTerm(u),
            front = deepcopy(mesh0),
            integrator = RK2(),
            redistributor = PoissonTangentialRedistributor(; iterations=2, pseudo_dt=0.08, omega=1.0, every=1),
        )
        integrate!(eq, 1.0; dt=dt)
        meshf = current_state(eq).mesh
        m = standard_curve_metrics(meshf; exact_area=exact_area, exact_perimeter=exact_perimeter)
        push!(rows, (
            case = "projected_rotation_circle",
            N = N,
            m...,
            E_kappa_inf = NaN,
            E_kappa_2 = NaN,
            E_n_bisector = NaN,
            E_n_osculating = NaN,
            E_n_laplacian = NaN,
            E_form = circle_shape_error(meshf; center=curve_vertex_centroid(meshf), R=R) / R,
        ))
    end
    return rows
end

function write_rows_csv(path::AbstractString, rows)
    isempty(rows) && return path
    names = propertynames(rows[1])
    open(path, "w") do io
        println(io, join(names, ","))
        for row in rows
            println(io, join((getproperty(row, n) for n in names), ","))
        end
    end
    return path
end

function print_rows(rows)
    for row in rows
        @printf("%-38s N=%4d E_A_arc=%.3e E_P_arc=%.3e E_k_inf=%.3e uniformity=%.3e E_form=%.3e\n",
                row.case, row.N, row.E_A_arc, row.E_P_arc, row.E_kappa_inf,
                row.uniformity, row.E_form)
    end
end

function run_intern_curve_benchmarks(;
    output = joinpath(@__DIR__, "output", "intern_curve_benchmarks.csv"),
    circle_Ns = (16, 32, 64, 128, 256),
    ellipse_Ns = (32, 64, 128),
    redistribution_Ns = (8, 16, 32, 64),
    projected_Ns = (8, 16, 32, 64),
)
    rows = NamedTuple[]
    append!(rows, run_static_circle_convergence(; Ns=circle_Ns))
    append!(rows, run_static_ellipse_convergence(; Ns=ellipse_Ns))
    append!(rows, run_redistribution_circle_study(; Ns=redistribution_Ns))
    append!(rows, run_projected_rotation_circle_study(; Ns=projected_Ns))

    isdir(dirname(output)) || mkpath(dirname(output))
    write_rows_csv(output, rows)
    print_rows(rows)
    println("Saved intern curve benchmark summary: $output")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_intern_curve_benchmarks()
end
