using Distributed
using PyPlot
using Parameters
using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Statistics
using NLsolve

# # Model Parameters
@with_kw mutable struct CRPar
    r = 1.0
    k = 3.0
    a = 1.0
    e = 0.7
    b = 0.8
    m = 0.4
    D = 0.0
    noise = 0.01
end

##Roz-Mac CR logistic Resource
function CR_mod!(du, u, p, t)
    @unpack r, k, a, b, e, m, D = p
    R, C, = u
    du[1] = r * R * (1 - R / k) - (a * R * C)/(R + b)
    du[2] = (e * a * R * C)/(R + b) - m * C - D * C^2
    return
end

##Stochastic CR Roz-Mac CR logistic Resource
function CR_mod(u, p)
    du = similar(u)
    CR_mod!(du, u, p, 0.0)
    return du
end

#Add noise to Stochastic
function stoch_CR_Mod!(du, u, p2, t)
    @unpack  noise = p2

    du[1] = noise * u[1]
    du[2] = noise * u[2]
    return du 
end

##Isocline plotting function
function retrieve_RC(vectordata, RC)
    newarray = zeros(size(vectordata)[1],size(vectordata)[2])

    for i in 1:size(vectordata)[1]
        for j in 1:size(vectordata)[2]
            newarray[i,j] = vectordata[i,j][RC]
        end
    end
    return newarray
end

##Plot Resource isocline
function res_isocline(R, p)
    @unpack r, k, a, b = p
    return (r/a) * (1 - R/k) * (R+b)
end

##Plot Consumer isocline
function cons_isocline(C, p)
    @unpack m, a, e, b, D = p
    denom = e * a - D * C - m
    if denom <= 0
        return NaN 
    else
        return (D * C * b + m) / denom
    end
end

##Plot different isoclines based on amount of dampening (D)
let
    par = CRPar()
    par1 = CRPar(D = 0.05)
    par2 = CRPar(D = 0.2)
    par3 = CRPar(D = 0.5)

    Rrange = 0.0:0.1:3.0
    Crange = 0.0:0.1:3.0
    resconrange = 0.0:0.005:3.0

    CR_mod_iso_plot = figure()
    plot(collect(Rrange), [res_isocline(R, par) for R in Rrange], color="navy")
    plot([cons_isocline(C, par) for C in Crange], collect(Crange),color="#FF9191")
    plot([cons_isocline(C, par1) for C in Crange], collect(Crange),color="#DC143C")
    plot([cons_isocline(C, par2) for C in Crange], collect(Crange),color="#B22222")
    plot([cons_isocline(C, par3) for C in Crange], collect(Crange),color="#800000")

    xlabel("Resource")
    ylabel("Consumer")
    xlim(0,3)
    ylim(0,3)

    #savefig("/Users/reillyoconnor/Desktop/Julia Projects/Slow-Fast/Julia/Figures/Isoclines.pdf", dpi = 300)
    return CR_mod_iso_plot
end

#Functions to calculate real and imaginary eigenvalues
# maximum(real.(eigvals(M))), where eigvals is from the standard library LinearAlgebra
λ_stability(M) = maximum(real.(eigvals(M)))
imag_eig(M) = maximum(imag(eigvals(M)))

calc_λ1(eq, par) = λ_stability(ForwardDiff.jacobian(eq -> CR_mod(eq, par), eq))

calc_λe(eq, par) = imag_eig(ForwardDiff.jacobian(eq -> CR_mod(eq, par), eq))

#Test that each component of local stability analysis works
par = CRPar()
tspan = (0.0, 1000.0)
u0 = [2.0, 1.0]

prob = ODEProblem(CR_mod!, u0, tspan, par)

sol = solve(prob)

eq = nlsolve((du, u) -> CR_mod!(du, u, par, 0.0), sol.u[end]).zero

###########################################################################
##### Investigate Aspects of stability over changing a (attack rate) ######
###########################################################################

##Real and Imaginary Eigenvalues
let
    avals = 0.60:0.01:1.3
    stab = fill(NaN, 2, length(avals))
    par = CRPar(D = 0.0)

    for (i, a) in enumerate(avals)
        par.a = a
        print(par)
        prob = ODEProblem(CR_mod!, u0, tspan, par)
        sol = solve(prob)
        eq = nlsolve((du, u) -> CR_mod!(du, u, par, 0.0), sol.u[end]).zero
        stab[1, i] = calc_λ1(eq, par)
        stab[2, i] = calc_λe(eq, par)
    end

    CR_Damp_eq = figure()
    subplot(211)
    plot(avals, stab[1, :])
    axhline(0, color = "black", linestyle = "--")
    ylabel(L"\lambda_1",fontsize=16,fontweight=:bold)

    subplot(212)
    plot(avals, stab[2, :])
    xlabel("Attack Rate (a)",fontsize=16,fontweight=:bold)
    ylabel(L"\lambda_ imag",fontsize=16,fontweight=:bold)
    
    tight_layout()

    #savefig("/Users/reillyoconnor/Desktop/Julia Projects/Slow-Fast/Julia/Figures/Eigenvalues_D0.pdf", dpi = 300)

    return CR_Damp_eq
end


##Function for Bifurcation Diagram
    function minmaxbifurc(tend)
    a_range = 0.1:0.001:3.0
    minC = zeros(length(a_range))
    maxC = zeros(length(a_range))

    u0 = [0.5, 0.5]
    t_span = (0.0, tend)
    trans = tend - 100.0
    remove_transient = trans:1.0:tend

    for (ai, a_vals) in enumerate(a_range)
        p = CRPar(a = a_vals, D = 0.00)
        print(p)
        prob = ODEProblem(CR_mod!, u0, t_span, p)
        sol = DifferentialEquations.solve(prob, reltol = 1e-8)
        solend = sol(remove_transient)
        solend[2, 1:end]

        minC[ai] = minimum(solend[2, :])
        maxC[ai] = maximum(solend[2, :])
    end

    return hcat(a_range, minC, maxC)
end

##Bifurcation Diagram
let
    data = minmaxbifurc(1000.0)
    minmaxplot = figure()
    scatter(data[:,1], data[:,2], label="C1 min", marker="o", s=1)
    scatter(data[:,1], data[:,3], label="C1 min", marker="o", s=1)
    xlabel("Attack Rate (a)")
    ylabel("Consumer Max/Min")

    #savefig("/Users/reillyoconnor/Desktop/Julia Projects/Slow-Fast/Figures/bifurcations_D005.pdf", dpi = 300)

    return minmaxplot
end

## CV Figure - Fig 2B in Manuscript - Loop over changing ea/m (via a)
## in p = CRPar() alter D to desired level of dampening
let 
    amaxs = 0.80:0.01:3.0
    maxhold = fill(0.0,length(amaxs),2)
    stdhold = fill(0.0,length(amaxs),1)
    meanhold = fill(0.0,length(amaxs),1)
    cvhold = fill(0.0,length(amaxs),1)
    taylorhold = fill(0.0,length(amaxs),1)
    
    ## Lets say we want the solutions at only certain time steps (just make sure it is inside of `t_span`! extrapolation is the road to sadness)
    ts = range(0, 1000, length = 1000)
    t_span = (0.0, 1000.0)

    for i = 1:length(amaxs)
        p = CRPar(a = amaxs[i], noise = 0.1, D = 0.0)
        print(p)
        u0 = [1.0, 0.5]
        prob_stoch = SDEProblem(CR_mod!, stoch_CR_Mod!, u0, t_span, p)
        sol_stoch = solve(prob_stoch, reltol = 1e-15)
        grid_sol = sol_stoch(ts)
        grid_sol.t
        grid_sol.u
        maxhold[i,1] = amaxs[i]
        maxhold[i,2] = maximum(grid_sol[2,900:1000])
        stdhold[i] = std(grid_sol[2,900:1000])
        meanhold[i] = mean(grid_sol[2,900:1000])
        cvhold[i] = stdhold[i]/meanhold[i]
        taylorhold[i]= (stdhold[i]^2)/meanhold[i]
    end

    CV_figure = figure(figsize = (6,4))

    #subplot(311)
    #plot(maxhold[1:length(amaxs),1],stdhold[1:length(amaxs)], color="black", #linewidth=2)
    #ylabel("SD (C)",fontsize=12,fontweight=:bold)
    #ylim(0,0.6)


    #subplot(312)
    #plot(maxhold[1:length(amaxs),1],meanhold[1:length(amaxs)], color="black", linewidth=3)
    #ylabel("Mean (C)",fontsize=12,fontweight=:bold)
    #ylim(0, 1.4)

    #subplot(313)
    plot(maxhold[1:length(amaxs),1],cvhold[1:length(amaxs)], color="black")
    xlabel("Attack Rate (a)",fontsize=12,fontweight=:bold)
    ylabel("CV (C)",fontsize=12,fontweight=:bold)
    ylim(0, 2.25)
    xlim(0.8, 3.0)

    #subplot(414)
    #plot(maxhold[1:length(amaxs),1],taylorhold[1:length(amaxs)])
    #xlabel("Attack Rate (a)",fontsize=16,fontweight=:bold)
    #ylabel("Taylors Power Law (C)",fontsize=16,fontweight=:bold)

    tight_layout()

    #savefig("/Users/reillyoconnor/Desktop/Julia Projects/Slow-Fast/Julia/Figures/Fig 2B - CR CV (attack rate).jpeg", dpi = 300)

    return CV_figure
end
