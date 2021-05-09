#Ensemble of tracers

using .TracerAdvDiff_QG

using FourierFlows, Plots, Distributions
using FFTW: rfft, irfft
import GeophysicalFlows.MultiLayerQG

#Set up the MultiLayerQG.Problem to test with the modified module.

#Choose CPU or GPU
dev = CPU()

#Numerical and time-stepping parameters
nx = 64        # 2D resolution = nx^2
ny = nx

stepper = "FilteredRK4";  # timestepper
Δt = 0.01                 # timestep
nsubs  = 1                # number of time-steps for plotting (nsteps must be multiple of nsubs)
nsteps = 2000nsubs        # total number of time-steps


#Physical parameters for a two layer QG_problem
Lx = 2π        # domain size
μ = 5e-2       # bottom drag
β = 5          # the y-gradient of planetary PV

nlayers = 2     # number of layers
f0, g = 1, 1    # Coriolis parameter and gravitational constant
H = [0.2, 0.8]  # the rest depths of each layer
ρ = [4.0, 5.0]  # the density of each layer

U = zeros(nlayers) # the imposed mean zonal flow in each layer
U[1] = 1.0
U[2] = 0.0

#Set initial conditions.
ϵ = 0.3;
x, y = gridpoints(QG_prob.grid)

q_1_i = @.  ϵ * cos(4π / Lx * x_QG) * exp(-(x^2 + y^2) / 8)
q_2_i = @. -ϵ * cos(4π / Lx * x_QG) * exp(-(x^2 + y^2) / 8)

q_i = zeros(nx,ny,2)

q_i[:, :, 1] = q_1_i
q_i[:, :, 2] = q_2_i
qh_i = QG_prob.timestepper.filter .* rfft(q_i, (1, 2))         # only apply rfft in dims=1, 2
q_i  = irfft(qh_i, QG_prob.grid.nx, (1, 2))                    # only apply irfft in dims=1, 2

#=
    Now, instead of one `TracerAdvDiff_QG.Problem`, initialise an array of `TracerAdvDiff_QG.Problem`(s). At the moment each AD_prob has a corresponding
    QG_prob (mainly as a result of how I wrote the delay time function). There may be way to only use one flow but for now this is how I set it up.
=#
#Set diffusivity
κ = 0.01

#Set number of tracer problems, delay time of placing tracer into the flow, and blank arrays for shortcuts and clock time for QG whe tracer is dropped in
n = 10
AD_probs = Array{FourierFlows.Problem}(undef, n)
QG_probs = Array{FourierFlows.Problem}(undef, n)
delay_times = range(0, 10, step = 1)
sol_AD, cl_AD, v_AD, p_AD = Array{Any}(undef, n), Array{Any}(undef, n), Array{Any}(undef, n), Array{Any}(undef, n)

for i in 1:n
    #Set the tracer advection probelms with different delay times and the corresponding QG problems.
    QG_probs[i] = MultiLayerQG.Problem(nlayers, dev; nx=nx, Lx=Lx, f₀=f0, g=g, H=H, ρ=ρ, U=U, dt=Δt, stepper=stepper, μ=μ, β=β)
    MultiLayerQG.set_q!(QG_probs[i], q_i)
    AD_probs[i] = TracerAdvDiff_QG.Problem(;prob = QG_probs[i], delay_time = delay_times[i], nsubs = nsubs, kap = κ)
    sol_AD[i], cl_AD[i], v_AD[i], p_AD[i] = AD_probs[i].sol, AD_probs[i].clock, AD_probs[i].vars, AD_probs[i].params
end

#The grid is the same for all so just define one grid (can change this if need be).
g_AD = AD_probs[1].grid
x_AD, y_AD = gridpoints(g_AD)
x, y = g_AD.x, g_AD.y

#Set Gaussian blob initial condition in all problems
μIC = [0, 0]
Σ = [1 0; 0 1]
blob = MvNormal(μIC, Σ)
blob_IC(x, y) = pdf(blob, [x, y])
C₀ = @. blob_IC(x_AD, y_AD)

for i in 1:n
    TracerAdvDiff_QG.QGset_c!(AD_probs[i], C₀)
end

#Define blank arrays to save plots.
lower_layer_tracer_plots_AD = Plots.Plot{Plots.GRBackend}[]
upper_layer_tracer_plots_AD = Plots.Plot{Plots.GRBackend}[]

#Step the tracer advection problems forward and plot at the desired time step.
for i in 1:n
    #Define frequency at which to save a plot.
    #plot_time_AD is when to get the first plot, plot_time_inc is at what interval subsequent plots are created.
    #Setting them the same gives plots at equal time increments.
    plot_time_AD, plot_time_inc = 0.2, 0.2
    while cl_AD[i].step <= nsteps
        if cl_AD[i].step == 0
            tp_u = heatmap(x, y, v_AD[i].c[:, :, 1]',
                    aspectratio = 1,
                    c = :balance,
                    xlabel = "x",
                    ylabel = "y",
                    colorbar = true,
                    xlim = (-g_AD.Lx/2, g_AD.Lx/2),
                    ylim = (-g_AD.Ly/2, g_AD.Ly/2),
                    title = "C(x,y,t), t = "*string(round(cl_AD[i].t; digits = 2)));
            push!(upper_layer_tracer_plots_AD, tp_u)
            tp_l = heatmap(x, y, v_AD[i].c[:, :, 2]',
                    aspectratio = 1,
                    c = :balance,
                    xlabel = "x",
                    ylabel = "y",
                    colorbar = true,
                    xlim = (-g_AD.Lx/2, g_AD.Lx/2),
                    ylim = (-g_AD.Ly/2, g_AD.Ly/2),
                    title = "C(x,y,t), t = "*string(round(cl_AD[i].t; digits = 2)))
            push!(lower_layer_tracer_plots_AD, tp_l)
        elseif round(Int64, cl_AD[i].step) == round(Int64, plot_time_AD*nsteps)
            tp_u = heatmap(x, y, v_AD[i].c[:, :, 1]',
                    aspectratio = 1,
                    c = :balance,
                    xlabel = "x",
                    ylabel = "y",
                    colorbar = true,
                    xlim = (-g_AD.Lx/2, g_AD.Lx/2),
                    ylim = (-g_AD.Ly/2, g_AD.Ly/2),
                    title = "C(x,y,t), t = "*string(round(cl_AD[i].t; digits = 2)))
            push!(upper_layer_tracer_plots_AD, tp_u)
            tp_l = heatmap(x, y, v_AD[i].c[:, :, 2]',
                    aspectratio = 1,
                    c = :balance,
                    xlabel = "x",
                    ylabel = "y",
                    colorbar = true,
                    xlim = (-g_AD.Lx/2, g_AD.Lx/2),
                    ylim = (-g_AD.Ly/2, g_AD.Ly/2),
                    title = "C(x,y,t), t = "*string(round(cl_AD[i].t; digits = 2)))
            push!(lower_layer_tracer_plots_AD, tp_l)
            plot_time_AD += plot_time_inc
        end
        stepforward!(AD_probs[i], nsubs)
        TracerAdvDiff_QG.QGupdatevars!(AD_probs[i])
        #Updates the velocity field in advection problem to the velocity field in the MultiLayerQG.Problem at each timestep.
        TracerAdvDiff_QG.vel_field_update!(AD_probs[i], QG_probs[i], nsubs)
    end
end

plot_top = plot(upper_layer_tracer_plots_AD[1], upper_layer_tracer_plots_AD[2], 
                upper_layer_tracer_plots_AD[3], upper_layer_tracer_plots_AD[4],
                upper_layer_tracer_plots_AD[5], upper_layer_tracer_plots_AD[6])
     
#Display the tracer advection in the lower layer.
plot_bottom = plot(lower_layer_tracer_plots_AD[1], lower_layer_tracer_plots_AD[2], 
                   lower_layer_tracer_plots_AD[3], lower_layer_tracer_plots_AD[4],
                   lower_layer_tracer_plots_AD[5], lower_layer_tracer_plots_AD[6])

plot(plot_top, plot_bottom, layout=(2, 1), size=(1200, 1200))

plot_top = plot(upper_layer_tracer_plots_AD[7], upper_layer_tracer_plots_AD[8], 
                upper_layer_tracer_plots_AD[9], upper_layer_tracer_plots_AD[10],
                upper_layer_tracer_plots_AD[11], upper_layer_tracer_plots_AD[12])
     
#Display the tracer advection in the lower layer.
plot_bottom = plot(lower_layer_tracer_plots_AD[7], lower_layer_tracer_plots_AD[8], 
                   lower_layer_tracer_plots_AD[9], lower_layer_tracer_plots_AD[10],
                   lower_layer_tracer_plots_AD[11], lower_layer_tracer_plots_AD[12])

plot(plot_top, plot_bottom, layout=(2, 1), size=(1200, 1200))

#There are more plots though at the moment only the first two problems' plots are displayed