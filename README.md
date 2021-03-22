# QG_tracer_advection
## Passive tracer advection using quasigeostrophic flow

Hi all!

I have set up this repository so you can see the modifcations I have made to the tracer advection module.
It can now take in a MultiLayerQG.Problem and use the QG flow the problem produces to advect a passive tracer (at least that is what it looks like it is doing)!
This is done in TracerAdvDiff_QG.jl which is a copy of the module already available with my additions to make it work with a QG flow.
There is an example (test_TrAdvDiff_QG.jl) that shows how it works and produces plots of the tracer advection in the upper and lower layer.
If you download both and load the module TracerAdvDiff_QG.jl into the local environment it should all run fine.

Be great to have your thoughts!

If it looks to be doing the right thing I am happy to add it to the PassiveTracerAdvection module (provided Navid and Greg are happy about this!)

Thanks very much!



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
