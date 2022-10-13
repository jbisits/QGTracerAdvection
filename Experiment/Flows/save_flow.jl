include("NewParamsSquare/Square_new_params_256domain_512res.jl")

nsubs  = 1              #Set the number of steps the simulation takes at each iteration.      
nsteps = 18000          #Set the total amount of time steps the advection-diffusion simulation should run for.

#Set delay times (that is flow for some length of time, then drop tracer in)
delay_time = Δt̂ * 6000
#delay_time = 0
#Set the frequency at which to save data
save_freq = 100

total_steps = nsteps + 6000

function get_u(prob)
    sol, params, vars, grid = prob.sol, prob.params, prob.vars, prob.grid
  
    @. vars.qh = sol
    streamfunctionfrompv!(vars.ψh, vars.qh, params, grid)
    @. vars.uh = -im * grid.l * vars.ψh
    invtransform!(vars.u, vars.uh, params)
  
    return vars.u .+ params.U
end

filename = "zonal_velocty.jld2"
out = Output(QGProb, filename, (:u, get_u))

saveproblem(out)

#Simulation loop
while QGClock.step <= total_steps

    if QGClock.step % save_freq == 0
        saveoutput(out)

        log = @sprintf("step: %04d", QGClock.step)
        println(log)
    end
    stepforward!(QGProb, nsubs)
    MultiLayerQG.updatevars!(QGProb)

end

## Zonal velocity

using JLD2, CairoMakie

saved_u = jldopen(joinpath(pwd(), "Flows/zonal_velocty.jld2"))
keys(saved_u)

iterations = parse.(Int, keys(saved_u["snapshots/t"]))
t = [saved_u["snapshots/t/$i"] for i ∈ iterations]
# The flow was spun up for 6000 time steps then tracer released, so we only want from 6000 onwards
# to get the flow that advected the tracer.

find_flow = findall(iterations .>= 6000)
flow_iterations = iterations[find_flow]

u_ts_upper = [saved_u["snapshots/u/$i"][:, :, 1] for i ∈ flow_iterations]
u_ts_lower = [saved_u["snapshots/u/$i"][:, :, 2] for i ∈ flow_iterations]

x, y = saved_u["grid/x"], saved_u["grid/y"]

close(saved_u)

## Time mean zonal flow
time_mean_zonal_flow_upper = mean(u_ts_upper)
time_mean_zonal_flow_lower = mean(u_ts_lower)

fig1 = Figure()
titles = ["Time mean zonal velocity (upper layer)", "Time mean zonal velocity (lower layer)"]
ax = [Axis(fig1[1, i],
        title = titles[i],
        xlabel = "x̂",
        ylabel = "ŷ",
        aspect = 1) for i ∈ 1:2]
hm = CairoMakie.heatmap!(ax[1], x, y, time_mean_zonal_flow_upper)
Colorbar(fig1[2, 1], hm, vertical = false, label = "û", flipaxis = false)
hm = CairoMakie.heatmap!(ax[2], x, y, time_mean_zonal_flow_lower)
Colorbar(fig1[2, 2], hm, vertical = false, label = "û", flipaxis = false)
fig1
save("tmean_zonalvelocity.png", fig1)
## Time and zonal mean of the zonal flow

time_zonal_mean_upper = reshape(mean(time_mean_zonal_flow_upper; dims = 1), :)
time_zonal_mean_lower = reshape(mean(time_mean_zonal_flow_lower; dims = 1), :)

fig2 = Figure()
titles = ["Time and zonal mean\nzonal velocity (upper layer)", "Time and zonal mean\nzonal velocity (lower layer)"]
ax = [Axis(fig2[1, i],
        title = titles[i],
        xlabel = "ŷ",
        ylabel = "u̅",
        aspect = 1) for i ∈ 1:2]

lines!(ax[1], y, time_zonal_mean_upper)
lines!(ax[2], y, time_zonal_mean_lower)        
fig2
save("tzmean_zonalvelocity.png", fig2)