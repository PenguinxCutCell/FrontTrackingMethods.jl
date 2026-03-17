# Running tests

Run the package test-suite from the project root:

```bash
cd /path/to/FrontTrackingMethods.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

Notes and tips
- The test harness uses local `FrontIntrinsicOps.jl` (see `Project.toml` / `Manifest.toml`).
- If you see warnings about replaced modules while developing tests, ensure
  `test/test_utils.jl` is only included once (the test dispatcher `test/runtests.jl`
  should include `test_utils.jl` and call `using .TestUtils`).
- Long convergence/benchmark runs are available in `test/` and may take
  several minutes depending on resolutions and mesh refinements.
