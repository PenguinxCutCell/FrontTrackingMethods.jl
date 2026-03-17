module ExampleUtils

export ensure_dir, default_times, output_dir_for

function ensure_dir(path::AbstractString)
    isdir(path) || mkpath(path)
    return path
end

function output_dir_for(example_name::AbstractString)
    return ensure_dir(joinpath(@__DIR__, "output", example_name))
end

function default_times(tf::Real, nframes::Int=120)
    nframes ≥ 2 || error("nframes must be at least 2")
    return collect(range(0.0, stop=Float64(tf), length=nframes))
end

end
