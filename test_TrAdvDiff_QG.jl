#Test TracerAdvDiff_copy.jl using a two layer QG flow (from the geophysical flows package).

include("TracerAdvDiff_QG.jl")

using .TracerAdvDiff_QG

using FourierFlows, Plots
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

#Setup QG_problem and make easier to access the parts of the struct
QG_prob = MultiLayerQG.Problem(nlayers, dev; nx=nx, Lx=Lx, f₀=f0, g=g, H=H, ρ=ρ, U=U, dt=Δt, stepper=stepper, μ=μ, β=β)

sol_QG, cl_QG, pr_QG, vs_QG = QG_prob.sol, QG_prob.clock, QG_prob.params, QG_prob.vars
x_QG, y_QG = QG_prob.grid.x, QG_prob.grid.y

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

MultiLayerQG.set_q!(QG_prob, q_i)

#Set diffusivity
κ = 0.01
#Set the tracer advection probelm by passing in the QG problem 
AD_prob = TracerAdvDiff_QG.Problem(;prob = QG_prob, kap = κ)
sol_AD, cl_AD, v_AD, p_AD, g_AD = AD_prob.sol, AD_prob.clock, AD_prob.vars, AD_prob.params, AD_prob.grid
x_AD, y_AD = gridpoints(g_AD)
#Set the (same) initial condition in both layers.
a = 20
C₀_func(x, y) = log(1 + cosh(a)^2 / (cosh(a * sqrt(x^2 + y^2))^2)) / (2 * a)

C₀ = @. C₀_func(x_AD, y_AD)

TracerAdvDiff_QG.QGset_c!(AD_prob, C₀)

#Plot of initial condition in the upper layer.
heatmap(x_AD[:, 1], y_AD[1, :], v_AD.c[:, :, 1],
        title = "Initial tracer concentration",
        xlabel = "x",
        ylabel = "y",
        color = :balance,
        aspecetratio = 1)

#Define blank arrays in which to store the plots of tracer diffusion in each layer.
lower_layer_tracer_plots_AD = Plots.Plot{Plots.GRBackend}[]
upper_layer_tracer_plots_AD = Plots.Plot{Plots.GRBackend}[]
#Define frequency at which to save a plot.
#plot_time_AD is when to get the first plot, plot_time_inc is at what interval subsequent plots are created.
#Setting them the same gives plots at equal time increments. (Might be a better work around)
plot_time_AD, plot_time_inc = 0.2, 0.2
#Step the tracer advection problem forward and plot at the desired time step.

aspectratio = 1,
c = :balance,
xlabel = "x",
ylabel = "y",

kwargs = (
         xlabel = "x",
         ylabel = "y",
    aspectratio = 1,
          color = :balance,
       colorbar = true,
           xlim = (-g_AD.Lx/2, g_AD.Lx/2),
           ylim = (-g_AD.Ly/2, g_AD.Ly/2)
)


while cl_AD.step <= nsteps
    if cl_AD.step == 0
        tp_u = heatmap(x_AD[:, 1], y_AD[1, :], v_AD.c[:, :, 1],
        aspectratio = 1,
        c = :balance,
        xlabel = "x",
        ylabel = "y",
        title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)));
        push!(upper_layer_tracer_plots_AD, tp_u)
        tp_l = heatmap(x_AD[:, 1], y_AD[1, :], v_AD.c[:, :, 2],
            aspectratio = 1,
            c = :balance,
            xlabel = "x",
            ylabel = "y",
            title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)))
        push!(lower_layer_tracer_plots_AD, tp_l)
    elseif round(Int64, cl_AD.step) == round(Int64, plot_time_AD*nsteps)
        tp_u = heatmap(x_AD[:, 1], y_AD[1, :], v_AD.c[:, :, 1],
        aspectratio = 1,
        c = :balance,
        xlabel = "x",
        ylabel = "y",
        title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)))
        push!(upper_layer_tracer_plots_AD, tp_u)
        tp_l = heatmap(x_AD[:, 1], y_AD[1, :], v_AD.c[:, :, 2],
            aspectratio = 1,
            c = :balance,
            xlabel = "x",
            ylabel = "y",
            title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)))
        push!(lower_layer_tracer_plots_AD, tp_l)
        global plot_time_AD += plot_time_inc
    end
    stepforward!(AD_prob, nsubs)
    TracerAdvDiff_QG.QGupdatevars!(AD_prob)
    #Updates the velocity field in advection problem to the velocity field in the MultiLayerQG.Problem at each timestep.
    vel_field_update!(AD_prob, QG_prob, nsubs)
end
#Need to set this up so this does not need to be hardcoded.
#Display the tracer advection in the upper layer.
plot_top = plot(upper_layer_tracer_plots_AD[1], upper_layer_tracer_plots_AD[2], 
                upper_layer_tracer_plots_AD[3], upper_layer_tracer_plots_AD[4],
                upper_layer_tracer_plots_AD[5], upper_layer_tracer_plots_AD[6])
     
#Display the tracer advection in the lower layer.
plot_bottom = plot(lower_layer_tracer_plots_AD[1], lower_layer_tracer_plots_AD[2], 
                   lower_layer_tracer_plots_AD[3], lower_layer_tracer_plots_AD[4],
                   lower_layer_tracer_plots_AD[5], lower_layer_tracer_plots_AD[6])

plot(plot_top, plot_bottom, layout=(2, 1), size=(1200, 1200))
