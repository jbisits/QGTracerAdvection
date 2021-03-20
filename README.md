# QG_tracer_advection
## Passive tracer advection using quasigeostrophic flow

Hi all!

I have set up this repository so you can see the modifcations I have made to the tracer advection module.
It can now take in a MultiLayerQG.Problem and use the QG flow the problem produces to advect a passive tracer (at least that is what it looks like it is doing)!
This is done in TracerAdvDiff_copy.jl which is a copy of the module already available with my additions to make it work with a QG flow.
There is an example (test_trcopy.jl) that shows how it works and produces plots of the tracer advection in the upper and lower layer.
If you download both and load the module TracerAdvDiff_copy.jl into the local environment it should all run fine.

Be great to have your thoughts!

If it looks to be doing the right thing I am happy to add it to the PassiveTracerAdvection module (provided Navid and Greg are happy about this!)

Thanks very much!
