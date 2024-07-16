using Distributions
using Parameters
using Statistics
using PyPlot
using DataFrames
using CSV

# # Figure 1 a
## Delayed ricker logistic population model for growth rate changing for
## Deterministic (eigenvalue and min/max), and stochastic (CV)

## Set Parameters for Model
@with_kw mutable struct RickerPar
    r = 1.0
    K = 1.0
    σ = 0.1
end

##Eigenvalue stability
λ1_stability(r) = 1 - r

##Deterministic Ricker Model
function ricker(N0, n_steps, p)
    @unpack r, K = p
    N = fill(0.0, n_steps)
    N[1] = N0
    for t in 1:(length(N) - 1)
        N[t + 1] = N[t] * exp(r * N[t] * (1 - N[t] / K))
    end
    return N
end

##Stochastic Ricker Model
function ricker_stochastic(N0, n_steps, p)
    @unpack r, K = p
    N = fill(0.0, n_steps)
    N[1] = N0
    for t in 1:(length(N) - 1)
        N[t + 1] = max(0.0, N[t] * exp(r * N[t] * (1 - N[t] / K)) + rand(Normal(0, p.σ)))
    end
    return N
end

## Bifuraction Function
function min_max(ts; tol = 1e-6)
    mins = eltype(ts)[]
    maxs = eltype(ts)[]
    for i in 2:(length(ts) - 1)
        if isapprox(ts[i - 1], ts[i], atol = tol)
            # equilbrium add to both
            if all(isapprox.(mins, ts[i], atol = tol))
                push!(mins, ts[i])
            end
            if all(isapprox.(maxs, ts[i], atol = tol))
                push!(maxs, ts[i])
            end
        elseif ts[i - 1] > ts[i] && ts[i] < ts[i + 1]
            # don't add nearby solutions, this is stupid ineffecient
            if all(isapprox.(mins, ts[i], atol = tol))
                push!(mins, ts[i])
            end
        elseif ts[i - 1] < ts[i] && ts[i] > ts[i + 1]
            # don't add nearby solutions, this is stupid ineffecient
            if all(isapprox.(maxs, ts[i], atol = tol))
                push!(maxs, ts[i])
            end
        end
    end
    return (mins = mins, maxs = maxs)
end

#Calculate CV Function
function CV(p::RickerPar)
    N0 = 0.3
    t_end = 10000
    t = 1:t_end
    t_attract = round(Int, length(t) / 2)
    N = ricker_stochastic(N0, t_end, p)
    return std(N[t_attract:end]) / mean(N[t_attract:end])
end

#Evaluate CV over changing R function
function CV(r_vals)
    p = RickerPar()
    CVs = fill(0.0, length(r_vals))
    for (i, r) in enumerate(r_vals)
        p.r = r
        CVs[i] = CV(p)
    end
    return CVs
end

#Plot Bifuraction Function
function plot_min_max()
    par = RickerPar()
    N0 = 0.3

    t_end = 5000
    r_vals = range(0.1, 3.1, length = 10000)

    for r in r_vals
        par.r = r
        N = ricker(N0, t_end, par)
        min_maxs = min_max(N[round(Int, t_end / 2):end], tol = 1e-3)

        plot(fill(r, length(min_maxs.mins)), min_maxs.mins, "k,")
        plot(fill(r, length(min_maxs.maxs)), min_maxs.maxs, "k,")
    end

    return
end

## Bifurcation diagram over changing rmax - Fig. 2A
let
    r_vals = range(0.5, 2.8, length = 1000)

    figure()

    # Eigenvalues
    #subplot(311)
    #plot(r_vals, λ1_stability.(r_vals), "k")
    #xlabel("Growth rate (r)")
    #ylabel(L"\lambda_1")

    # Min / Max
    #subplot(312)
    plot_min_max()
    xlabel("Growth rate (r)")
    ylabel("Min/Max")
    xlim(0.5, 3.0)

    # CV
    #subplot(313)
    #plot(r_vals, CV(r_vals), "k")
    #xlabel("Growth rate (r)")
    #ylabel("CV (σ/μ)")
    #xlim(0.6, 3.0)

    #savefig("/Users/reillyoconnor/Desktop/Julia Projects/Slow-Fast/Julia/Figures/Fig 2A - Ricker Bifuraction.jpeg", dpi = 300)

    tight_layout()

end

PyPlot.display_figs()
