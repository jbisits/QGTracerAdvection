#Example of a time dependent flow advection using the module TracerAdvDiff_QG

#The TracerAdvDiff_QG.jl file must first be loaded into the local environment. See README for instructions.

using .TracerAdvDiff_QG
using Plots, FourierFlows, Distributions
#This can also be run using the package `PassiveTracerFlows` using the below command
#import PassiveTracerFlows.TracerAdvectionDiffusion

dev = CPU()

#Numerical and time-stepping parameters
nx = 64        # 2D resolution = nx^2
ny = nx
stepper = "FilteredRK4";  # timestepper
Δt = 0.01                 # timestep
nsubs  = 1                # number of time-steps for plotting (nsteps must be multiple of nsubs)
nsteps = 2000nsubs        # total number of time-steps

#Domain and diffusivity
Lx = 2π
nx = 64
κ = 0.01

#u and v velcocity functions which depend on time t
uvel(x, y, t) = -(t/2)*sin(x)*cos(y) + 0.5*cos(x)*sin(y)
vvel(x, y, t) =  (t/2)*cos(x)*sin(y) - 0.5*sin(y)*cos(x)

prob = TracerAdvDiff_QG.Problem(; steadyflow=false, nx=nx, Lx=Lx, κ=κ, u=uvel, v=vvel, dt=Δt, stepper=stepper)
sol, cl, v, p, g = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
x_grid, y_grid = gridpoints(g)
x, y = g.x, g.y
#A Gaussian "blob" initial condition centred at μIC 
μIC = [0, 0]
Σ = [0.5 0; 0 0.5]
blob = MvNormal(μIC, Σ)
blob_IC(x, y) = pdf(blob, [x, y])
C₀ = @. blob_IC(x_grid, y_grid)

TracerAdvDiff_QG.set_c!(prob, C₀)

#=
If using PassiveTracerFlows use the following problem set up.
prob = TracerAdvectionDiffusion.Problem(; steadyflow=true, nx=nx, Lx=Lx, κ=κ, u=uvel, v=vvel, dt=Δt, stepper=stepper);
sol, cl, v, p, g = prob.sol, prob.clock, prob.vars, prob.params, prob.grid;
TracerAdvectionDiffusion.set_c!(prob,C₀)
=#

#Step the problem forward and visualise the output.
tracer_plots = Plots.Plot{Plots.GRBackend}[]
plot_time = 0.2
while cl.step <= nsteps
    if cl.step == 0
        tp = heatmap(x, y, v.c',
            aspectratio = 1,
            c = :balance,
            xlabel = "x",
            ylabel = "y",
            colorbar = true,
            xlim = (-g.Lx/2, g.Lx/2),
            ylim = (-g.Ly/2, g.Ly/2),
            title = "C(x,y,t), t = "*string(round(cl.t; digits = 2)))
        push!(tracer_plots, tp);
    elseif round(Int64, cl.step) == round(Int64, plot_time*nsteps)
        tp = heatmap(x, y, v.c',
            aspectratio = 1,
            c = :balance,
            xlabel = "x",
            ylabel = "y",
            colorbar = true,
            xlim = (-g.Lx/2, g.Lx/2),
            ylim = (-g.Ly/2, g.Ly/2),
            title = "C(x,y,t), t = "*string(round(cl.t; digits = 2)))
        push!(tracer_plots, tp)
        global plot_time = 0.2 + plot_time
    end
    stepforward!(prob, nsubs);
    TracerAdvDiff_QG.updatevars!(prob)
    #If using PassiveTracerFlows comment out above line and use line below
    #TracerAdvectionDiffusion.updatevars!(prob)
end
plot(tracer_plots..., size=(1200, 1200))
