# QG tracer advection

## Passive tracer advection using quasigeostrophic flow

The module `TracerAdvDiff_QG` builds on the pre-existing module `TraverAdvDiff` from the [`PassiveTracerFlows`](https://fourierflows.github.io/PassiveTracerFlowsDocumentation/v0.1.0/)[1] package.
The `TracerAdvDiff_QG` module sets up a `TracerAdvDiff_QG.Problem` from a `MultiLayerQG.Problem` (from the [`GeophysicalFlows`](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/)[2] package) and advects a passive tracer in the quasigeosrophic flow the `MultiLayerQG.Problem` generates.
An initial condition for tracer concentration is set in all layers then either a series of plots or a movie (or both) are generated showing the tracer advection.

The tracer advection in either a steady flow or a time dependent flow remains unchanged from the original version of the module `TraverAdvDiff`.

## Examples

### Tracer advection with QG flow 
An example of how the module advects a passive tracer in a QG flow is given in the script `test_TracerAdvDiff.jl`. 
Here the `MultiLayerQG.Problem` from the [documentation](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/generated/multilayerqg_2layer/) is used as the flow to advect the tracer.
There are two different initial conditions in the script: a Gaussian "blob" or a Guassian "strip".
Both can be moved to different locations on the grid by altering the mean and the concentration about the mean is changed by altering the variance.
The same initial condition is set in both layers then the problem is stepped forward and the tracer is advected.
The output is then a series of plots showing the tracer advection.

### Tracer advection with a steady flow
An example of how the module advects a passive tracer in a steady flow is given in the script `steadyflow_example.jl`.
In this case only one `Problem` is initialised but the velocity fields `u` and `v` must be provided as functions of `x` and `y`.
To run this example the `TracerAdvDiff_QG` module must be first loaded into the local environment using
```julia
julia> include("TracerAdvDiff_QG.jl")
```
**Note this script can also be run using the `PassiveTracerFlows` package.**
There are instructions in the script as to how to do this.

### Tracer advection with a time dependent flow
An example of how the module advects a passive tracer in a time dependent flow is given in the script `timedep_example.jl`.
In this case only one `Problem` is initialised but the velocity fields `u` and `v` must be provided as functions of `x`, `y` and `t`.
To run this example the `TracerAdvDiff_QG` module must be first loaded into the local environment using
```julia
julia> include("TracerAdvDiff_QG.jl")
```
**Note this script can also be run using the `PassiveTracerFlows` package.**
There are instructions in the script as to how to do this.
## Using the module

To run first clone the repository, e.g.,

```
git clone https://github.com/jbisits/QG_tracer_advection.git
```

Then while in the repository's local directory open Julia, activate, and instantiate the project

```julia
julia>]
(@v1.5) pkg> activate .
(QG_tracer_advection) pkg> instantiate
```

Then a backspace will bring you back to normal Julia's REPL:
```julia
julia>
```

You can now run the `test_TrAdvDiff_QG.jl` via
```julia
julia> include("test_TrAdvDiff_QG.jl")
```

Alternatively, after you've instantiated the project, you can run the `test_TrAdvDiff_QG.jl` script straight from the terminal via

```
$ julia --project test_TrAdvDiff_QG.jl
```
# References
[1] Gregory L. Wagner & Navid C. Constantinou. (2018). FourierFlows/FourierFlows.jl. Zenodo. https://doi.org/10.5281/zenodo.1161724

[2] Constantinou et al., (2021). GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs. _Journal of Open Source Software_, **6(60)**, 3053, doi:[10.21105/joss.03053](https://doi.org/10.21105/joss.03053)
