# QG tracer advection

This repository contains all the code to generate the tracer release simulations and analysis of the resulting data for the paper *"Can we accurately quantify a lateral diffusivity from a single tracer release?"* submitted to JPO in July 2022.

## NOTE: This repository has been archived as it is not actively maintained any more. Passive tracer advection-diffusion using a turbulent flow has now been implemented in [`PassiveTracerFlows.jl`](https://github.com/FourierFlows/PassiveTracerFlows.jl) as of version 0.6.1. `PassiveTracerFlows.jl` will be maintained so I recommend you use it instead of this code repository which I am not maintaining. 

## Passive tracer advection using quasigeostrophic flow

The module `TracerAdvDiff_QG` builds on the pre-existing module `TraverAdvDiff` from the [`PassiveTracerFlows.jl`](https://fourierflows.github.io/PassiveTracerFlowsDocumentation/stable/)(1) package.
The `TracerAdvDiff_QG` module sets up a `TracerAdvDiff_QG.Problem` from a `MultiLayerQG.Problem` (from the [`GeophysicalFlows`](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/)[2] package) and advects a passive tracer in the quasigeosrophic flow the `MultiLayerQG.Problem` generates.
An initial condition for tracer concentration is set in all layers then either a series of plots or a movie (or both) are generated showing the tracer advection.

The tracer advection in either a steady flow or a time dependent flow remains unchanged from the original version of the module `TraverAdvDiff`.

The module `MeasureMixing.jl` has functions that can be used to create plots or animations from saved simulation data, as well as various concentration diagnostics.

## Using the modules

To run first clone the repository, e.g.,

```julia
git clone https://github.com/jbisits/QG_tracer_advection.git
```

Then while in the repository's local directory open Julia, activate, and instantiate the project

```julia
julia>]
(@v1.6) pkg> activate .
(QG_tracer_advection) pkg> instantiate
```

Then a backspace will bring you back to normal Julia's REPL:

```julia
julia>
```

Before running the examples or using the modules they must first be loaded into the local environment using

```julia
julia> include("Modules/TracerAdvDiff_QG.jl")
```

and

```julia
julia> include("Modules/MeasureMixing.jl")
```

Alternatively, after you've instantiated the project and loaded the modules in, you can run the `QGflow_example.jl` script straight from the terminal via

```julia
julia --project Examples/QGflow_example.jl
```

## Examples

To run the examples the modules must first be loaded into the local envirnoment using the commands

```julia
julia> include("Modules/TracerAdvDiff_QG.jl")
```

and

```julia
julia> include("Modules/MeasureMixing.jl")
```

### Tracer advection with QG flow

An example of how the module advects a passive tracer in a QG flow is given in the script `QGflow_example.jl`.
Here the `MultiLayerQG.Problem` from the [GeophysicalFlows.jl documentation](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/generated/multilayerqg_2layer/) is used as the flow to advect the tracer.
There are two different initial conditions in the script: a Gaussian "blob" or a Gaussian "strip".
Both can be moved to different locations on the grid by altering the mean and the concentration about the mean is changed by altering the variance.
The same initial condition is set in both layers then the problem is stepped forward and the tracer is advected-diffused.
The concentration field is saved to a `.jld2` file at regular intervals specified by the user.
After the saved simulation data has been loaded it can then be passed to the various functions of the `MeasureMixing.jl` module as is shown in the example.

After the modules have been loaded into the local environment use the command

```julia
julia> include("Examples/QGflow_example.jl")
```

to run the example.

### Tracer advection with a steady flow

An example of how the module advects a passive tracer in a steady flow is given in the script `steadyflow_example.jl`.
In this case only one `Problem` is initialised but the velocity fields `u` and `v` must be provided as functions of `x` and `y`.

After the modules have been loaded into the local environment use the command

```julia
julia> include("Examples/steadyflow_example.jl")
```

to run the example.

**Note this script can also be run using the `PassiveTracerFlows` package.**
There are instructions in the script as to how to do this.

### Tracer advection with a time dependent flow

An example of how the module advects a passive tracer in a time dependent flow is given in the script `timedep_example.jl`.
In this case only one `Problem` is initialised but the velocity fields `u` and `v` must be provided as functions of `x`, `y` and `t`.

After the modules have been loaded into the local environment use the command

```julia
julia> include("Examples/timedep_example.jl")
```

to run the example.

**Note this script can also be run using the `PassiveTracerFlows` package.**
There are instructions in the script as to how to do this.

## Honours thesis

This folder contains the code that I will be using for my honours thesis that I (@jbisits) am currently undertaking at UNSW.

## References

[1] Constantinou, N. C. and Wagner, G. L. (2021). FourierFlows/PassiveTracerFlows.jl: PassiveTracerFlows v0.5.0 (Version v0.5.0). Zenodo. <https://doi.org/10.5281/zenodo.2535983>

[2] Constantinou et al., (2021). GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs. _Journal of Open Source Software_, **6(60)**, 3053, doi:[10.21105/joss.03053](https://doi.org/10.21105/joss.03053)
