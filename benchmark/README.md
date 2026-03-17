# Advection Benchmarks and Error Definitions

This folder contains benchmark drivers for interface-advection tests and computes:

- area error
- shape error (for circular interfaces)
- symmetric-difference area

The advection input used in these scripts is built from:

1. analytical velocity sampled on staggered faces of a uniform Cartesian mesh,
2. bilinear interpolation of those face values back to marker positions,

then passed into `AdvectionTerm`.

## Error formulas

For initial domain $A(0)$ and final domain $A(T)$:

- $E_{area} = \frac{|A(T)-A(0)|}{A(0)}$

For circular interfaces with center $(x_c,y_c)$ and radius $R$:

- $E_{shape}=\max_i |\mathrm{dist}(x_i)|$
- $\mathrm{dist}(x_i)=\sqrt{(x_i-x_c)^2+(y_i-y_c)^2}-R$

For two domains $A$ (initial) and $B$ (final):

- $A\triangle B = (A\cup B)\setminus (A\cap B)$
- $E_{sym}=|A\triangle B| = |A| + |B| - 2|A\cap B|$

`E_sym` is evaluated numerically on a uniform Cartesian quadrature grid over [0,1]².

## Cases

- `translation_circle.jl`
  - circle: $R=0.15$, center $(0.25,0.75)$
  - velocity: $(u,v)$ with $u=-v$ and piecewise reversal at $t=0.5T$
  - target: return to initial shape at $T=1$
  - timestep by $\mathrm{CFL}=u\,\Delta t/h=0.125$

- `solid_body_rotation_circle.jl`
  - circle: $R=0.15$, center $(0.5,0.75)$
  - velocity: $(u,v)=(2\pi(0.5-y), 2\pi(x-0.5))$
  - one period: $T=1$
  - timestep by $\mathrm{CFL}=u_{max}\,\Delta t/h=\pi/16\approx0.2$

- `zalesak_disk_rotation.jl`
  - notched disk: $R=0.15$, center $(0.5,0.75)$
  - notch width $0.05$, notch depth $0.25$
  - same velocity/timestep strategy as solid-body rotation
  - reported: $E_{area}$ and $E_{sym}$

- `single_vortex_circle.jl`
  - circle: $R=0.15$, center $(0.5,0.75)$
  - Rider–Kothe reversible vortex with
    $\phi=\pi^{-1}\sin^2(\pi x)\sin^2(\pi y)\cos(\pi t/T)$,
    $(u,v)=(\partial\phi/\partial y, -\partial\phi/\partial x)$
  - timestep by $\mathrm{CFL}=u_{max}\,\Delta t/h=0.125$ at $t=0$

## Run

From package root:

`julia --project=. benchmark/run_all.jl`

Results are saved to:

`benchmark/output/advection_errors.csv`

Run an individual case (with plots and animation):

`julia --project=. benchmark/single_vortex_circle.jl`

Per-case visuals are written in:

`benchmark/output/<case_name>/initial.png`, `final.png`, `animation.mp4`

## Marker-count convergence study

`convergence_markers.jl` sweeps each benchmark case over a range of marker
counts $N$ (default: 32, 64, 128, 256, 512), collects $E_{area}$, $E_{shape}$,
$E_{sym}$ vs.\ the mean edge length $h \approx 2\pi R / N$, estimates
convergence orders via a least-squares log-log fit, and saves results.

Run the full sweep (produces CSV + convergence PNGs):

`julia --project=. benchmark/convergence_markers.jl`

Or call from a Julia session with custom resolutions:

```julia
include("benchmark/convergence_markers.jl")
run_marker_convergence(resolutions=[32, 64, 128, 256], save_visuals=true)
```

Outputs are saved in:

- `benchmark/output/convergence/marker_convergence.csv`
- `benchmark/output/convergence/<case>_convergence.png`  (log-log plots)

### Convergence order estimator

Given pairs $(h_i, E_i)$, the order $p$ is estimated as the least-squares slope in log-log space:

$$
p = \frac{\sum_i (\ln h_i - \overline{\ln h})(\ln E_i - \overline{\ln E})}{\sum_i (\ln h_i - \overline{\ln h})^2}
$$

## EBIT conversion

Number of markers is defined by the number of crossing point between the front and a uniform Cartesian grid of spacing $h=1/N$.

For a **32×32 Cartesian mesh** on the unit square, the grid lines are at

$$
x=\frac{i}{32},\quad y=\frac{j}{32}, \qquad i,j=0,\dots,32.
$$

The circle is

$$
(x-0.25)^2+(y-0.75)^2=0.15^2.
$$

A marker is created each time the circle crosses a **mesh edge**, so we count intersections with:

* vertical grid lines (x=i/32),
* horizontal grid lines (y=j/32).

The circle spans:

* in (x): ([0.25-0.15,;0.25+0.15]=[0.10,;0.40]),
* in (y): ([0.75-0.15,;0.75+0.15]=[0.60,;0.90]).

Grid lines inside those intervals are:

* vertical: (x=\frac{4}{32},\frac{5}{32},\dots,\frac{12}{32}) → **9 lines**
* horizontal: (y=\frac{20}{32},\frac{21}{32},\dots,\frac{28}{32}) → **9 lines**

None of these are tangent, so each such line cuts the circle in **2 points**.

So:

$$
9\times 2 + 9\times 2 = 18+18 = 36.
$$

For a $2^k\times 2^k$ Cartesian mesh, the number of markers is:
- count vertical lines crossing the circle, multiplying by 2 for the two intersections per line,
- count horizontal lines crossing the circle, multiplying by 2 for the two intersections per line,
- sum these two contributions.

$$
M(N)=2\left(
\left\lceil \frac{2N}{5}\right\rceil-\left\lfloor \frac{N}{10}\right\rfloor-1
+
\left\lceil \frac{9N}{10}\right\rceil-\left\lfloor \frac{3N}{5}\right\rfloor-1
\right)
$$

| $k$ | Markers |
|----|---------|
| 2 | 4 |
| 3 | 12 |
| 4 | 20 |
| 5 | 36 |
| 6 | 76 |
| 7 | 156 |
| 8 | 308 |
| 9 | 612 |
| 10 | 1228 |
| 11 | 2460 |


