using DelimitedFiles

include("translation_circle.jl")
include("solid_body_rotation_circle.jl")
include("zalesak_disk_rotation.jl")
include("single_vortex_circle.jl")

function run_all_benchmarks(; resolutions=(128, 256))
    rows = Any[]

    for N in resolutions
        r = run_translation_circle(N=N, save_visuals=false)
        push!(rows, ("translation_circle", N, r.h, r.dt, r.CFL, r.E_area, r.E_shape, r.E_sym))
    end

    for N in resolutions
        r = run_solid_body_rotation_circle(N=N, save_visuals=false)
        push!(rows, ("solid_body_rotation_circle", N, r.h, r.dt, r.CFL, r.E_area, r.E_shape, r.E_sym))
    end

    for N_arc in resolutions
        r = run_zalesak_disk_rotation(N_arc=N_arc, save_visuals=false)
        push!(rows, ("zalesak_disk_rotation", N_arc, r.h, r.dt, r.CFL, r.E_area, NaN, r.E_sym))
    end

    for N in resolutions
        r = run_single_vortex_circle(N=N, save_visuals=false)
        push!(rows, ("single_vortex_circle", N, r.h, r.dt, r.CFL, r.E_area, r.E_shape, r.E_sym))
    end

    outdir = joinpath(@__DIR__, "output")
    isdir(outdir) || mkpath(outdir)
    outfile = joinpath(outdir, "advection_errors.csv")

    open(outfile, "w") do io
        println(io, "case,resolution,h,dt,cfl,E_area,E_shape,E_sym")
        for r in rows
            println(io, join(r, ","))
        end
    end

    println("Saved benchmark summary: $outfile")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_all_benchmarks()
end
