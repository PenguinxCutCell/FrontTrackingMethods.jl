# Installation

Install with the current Julia project activated. If you are developing
locally, add the package by path so local changes are used:

```julia
using Pkg
Pkg.develop(path=".")   # from the package root
Pkg.resolve()
```

To run the test-suite:

```bash
cd /path/to/FrontTrackingMethods.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

Notes:
- This package depends on `FrontIntrinsicOps.jl`. When running locally,
  ensure `FrontIntrinsicOps.jl` is available (the test setup uses a
  local path dependency in the Project/Manifest).
- The project uses `StaticArrays` and standard library modules.
