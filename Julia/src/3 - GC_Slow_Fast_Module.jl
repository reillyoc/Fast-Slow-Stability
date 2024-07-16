using Parameters
using ForwardDiff
using LinearAlgebra
using PyPlot
using DifferentialEquations
using NLsolve
using Statistics
using RecursiveArrayTools
using Noise
using Distributed
using StatsBase
using Random
using CSV
using DataFrames

## Here we will illustrate the case of a generalist predator that feeds on fast and slow channel
## Parameters are categorized by fast and slow macrohabitats/energy channels
## Parameters with "_2" indicate slow channel and those with "_1" indicaste fast channel

## Set parameters
@with_kw mutable struct GenMod_Par
    α_1 = 0.0   ##competitive influence of R1 on R2 - Set to 0
    α_2 = 0.0   ##competitve influence of R2 on R1 - Set to 0
    k_2 = 1.0
    k_1 = 1.0
    e_CR = 0.8
    e_PC = 0.6
    e_PR = 0.8
    m_P = 0.2

    a_PC_2 = 0.1
    a_PR_1 = 0.00 ## No omnivory - Set to 0
    a_PC_1 = 0.3
    a_PR_2 = 0.00 ## No omnivory - Set to 0
    h_PC = 0.60
    h_PR = 1.0
  
#slow CR
    r_2 = 1.0
    a_CR_2 = 0.5
    m_Cl = 0.10
    h_CRl = 1.0

#fast CR
    r_1 = 1.0
    a_CR_1 = 0.9
    m_Cp = 0.20
    h_CRp = 1.0

# noise
    noise = 0.01
end


## Generalist Module, with omnivory and resource competition
## Omnivory and competition set to 0 and thus result in equations given in manuscript - See Materials and Methods - eq, 6,7,8

function GenMod_model!(du, u, p, t)
    @unpack r_2, r_1, k_2, k_1, α_1, α_2, e_CR, e_PC, e_PR, a_CR_2, a_CR_1, a_PR_2, a_PR_1, h_CRl, h_CRp,h_PC, h_PR, m_Cl, m_Cp,m_P, a_PC_2,a_PC_1 = p 
    
    
    R_2, R_1, C_2, C_1, P = u
    
    du[1]= r_2 * R_2 * (1 - (α_1 * R_1 + R_2)/k_2) - (a_CR_2 * R_2 * C_2)/(1 + a_CR_2 * h_CRl * R_2) - (a_PR_2 * R_2 * P)/(1 + a_PR_2 * h_PR * R_2 + a_PR_1 * h_PR * R_1 + a_PC_2 * h_PC * C_2 + a_PC_1 * h_PC * C_1)
    
    du[2] = r_1 * R_1 * (1 - (α_2 * R_2 + R_1)/k_1) - (a_CR_1 * R_1 * C_1)/(1 + a_CR_1 * h_CRp * R_1) - (a_PR_1 * R_1 * P)/(1 + a_PR_2 * h_PR * R_2 + a_PR_1 * h_PR * R_1 + a_PC_2 * h_PC * C_2 + a_PC_1 * h_PC * C_1)

    du[3] = (e_CR * a_CR_2 * R_2 * C_2)/(1 + a_CR_2 * h_CRl * R_2) - (a_PC_2 * C_2 * P)/(1 + a_PR_2 * h_PR * R_2 + a_PR_1 * h_PR * R_1 + a_PC_2 * h_PC * C_2 + a_PC_1 * h_PC * C_1) - m_Cl * C_2
    
    du[4] = (e_CR * a_CR_1 * R_1 * C_1)/(1 + a_PC_1 * h_CRp * R_1) - (a_PC_1 * C_1 * P)/(1 + a_PR_2 * h_PR * R_2 + a_PR_1 * h_PR * R_1 + a_PC_2 * h_PC * C_2 + a_PC_1 * h_PC * C_1) - m_Cp * C_1

    du[5] = (e_PR * a_PR_2 * R_2 * P + e_PR * a_PR_1 * R_1 * P + e_PC * a_PC_2 * C_2 * P + e_PC * a_PC_1 * C_1 * P)/(1 + a_PR_2 * h_PR * R_2 + a_PR_1 * h_PR * R_1 + a_PC_2 * h_PC * C_2 + a_PC_1 * h_PC * C_1) - m_P * P

    return 
end


function GenMod_model(u, p)
    du = similar(u)
    GenMod_model!(du, u, p, 0.0)
    return du
end

## Adding stochasticity to model using gaussian white noise (SDEproblem)
function GenMod_stochmodel!(du, u, p2, t)
    @unpack  noise = p2

    du[1] = noise * u[1]
    du[2] = noise * u[2]
    du[3] = noise * u[3]
    du[4] = noise * u[4]
    du[5] = noise * u[5]

    return du 
end


#Functions to calculate real and imaginary eigenvalues
# maximum(real.(eigvals(M))), where eigvals is from the standard library LinearAlgebra
λ_stability(M) = maximum(real.(eigvals(M)))
imag_eig(M) = maximum(imag(eigvals(M)))

calc_λ1(eq, par) = λ_stability(ForwardDiff.jacobian(eq -> GenMod_model(eq, par), eq))

calc_λe(eq, par) = imag_eig(ForwardDiff.jacobian(eq -> GenMod_model(eq, par), eq))

#Test that each component of local stability analysis works
par = GenMod_Par()
tspan = (0.0, 1000.0)
u0 = [1.0, 1.0, 0.5, 0.5, 0.1]

prob = ODEProblem(GenMod_model!, u0, tspan, par)

sol = solve(prob)

eq = nlsolve((du, u) -> GenMod_model!(du, u, par, 0.0), sol.u[end]).zero

#stochastic C-R time series -- just to look at timeseries specifically whenever we want 
let
    u0 = [1.0, 1.0, 0.5, 0.5, 0.25]
    t_span = (0.0, 2000.0)
    rmax_add = 0.3
    par = GenMod_Par(noise = 0.01, α_1 = 0.0, α_2 = 0.0)

    par.a_CR_1 = (par.a_CR_1 + rmax_add)
    par.a_CR_2 = (par.a_CR_2 + rmax_add)
    par.a_PC_1 = (par.a_PC_1 + rmax_add)
    par.a_PC_2 = (par.a_PC_2 + rmax_add)

    prob_stoch = SDEProblem(GenMod_model!, GenMod_stochmodel!, u0, t_span, par)
    sol_stoch = solve(prob_stoch, reltol = 1e-15)
    ts_genmod = figure()
    plot(sol_stoch.t[1:end], sol_stoch[1, 1:end], label = "R1")
    plot(sol_stoch.t[1:end], sol_stoch[2, 1:end], label = "R2")
    plot(sol_stoch.t[1:end], sol_stoch[3, 1:end], label = "C2")
    plot(sol_stoch.t[1:end], sol_stoch[4, 1:end], label = "C1")
    plot(sol_stoch.t[1:end], sol_stoch[5, 1:end], label = "P")
    xlabel("Time")
    ylabel("Abundance")
    legend()
    return ts_genmod
end

/


###################################################################################
##### Investigate Aspects of stability over changing P & C1, C2 attack rates ######
###################################################################################

##### Eigenvalue Figure for increasing P & C1,C2 attack rates - Fig. 2C #####
## Deterministic Model
let
    rmax_increases = 0.01:0.01:0.60
    stab = fill(NaN, 2, length(rmax_increases))

    #Initialize your parameters
    par = GenMod_Par(noise = 0.00) 

    for (i, rmax_increase) in enumerate(rmax_increases)
        #Modify all attack rates by the percentage increase
        par.a_CR_1 = (par.a_CR_1 + rmax_increase)
        par.a_CR_2 = (par.a_CR_2 + rmax_increase)
        par.a_PC_1 = (par.a_PC_1 + rmax_increase)
        par.a_PC_2 = (par.a_PC_2 + rmax_increase)

        #Solve the ODE problem
        prob = SDEProblem(GenMod_model!, GenMod_stochmodel!, u0, tspan, par)
        sol = solve(prob, reltol = 1e-15)

        #Find the equilibrium state
        eq = nlsolve((du, u) -> GenMod_model!(du, u, par, 0.0), sol.u[end]).zero

        #Calculate and store the eigenvalues
        stab[1, i] = calc_λ1(eq, par)
        stab[2, i] = calc_λe(eq, par)

        #Reset the parameters to base values for next iteration
        par = GenMod_Par(noise = 0.00)
    end

    GenMod_eq = figure(figsize = (6,4))
    #subplot(211)
    plot(rmax_increases, stab[1, :], color = "black")
    axhline(0, color = "black", linestyle = "--")
    ylabel(L"\lambda_max", fontsize=16, fontweight=:bold)
    ylim(-0.175, 0.025)
    #axvspan(-0.40, 0.11, facecolor="gray", alpha=0.5)

    #subplot(212)
    #plot(rmax_increases, stab[2, :])
    #xlabel("Attack Rate (a)",fontsize=16,fontweight=:bold)
    #ylabel(L"\lambda_ imag",fontsize=16,fontweight=:bold)
    #axvspan(-0.40, 0.11, facecolor="gray", alpha=0.5)
    
    tight_layout()

    #savefig("/Users/reillyoconnor/Desktop/Julia Projects/Slow-Fast/Julia/Figures/Fig. 2C - Generalist Module Eigenvalue Checkmark.jpeg", dpi = 300)

    return GenMod_eq
end


##### Consquences of Accelerating Ecosystems (Fig. 4C .csv - see R code for fig.) #####
## Stochastic Model
let
    percent_increases = 0.0:0.001:0.05 #Increase in rmax via a - 0-5%
    rmax_high = 0.11
    rmax_low = -0.40
    n_seeds = 100
    maxhold = fill(0.0, length(percent_increases), 1)
    stdhold = fill(0.0, length(percent_increases), 2)
    meanhold = fill(0.0, length(percent_increases), 2)
    cv_high_all = fill(0.0, length(percent_increases), n_seeds)
    cv_low_all = fill(0.0, length(percent_increases), n_seeds)
    
    
    ts = range(0, 2000, length = 2000)
    t_span = (0.0, 2000.0)

    #Initialize your parameters
    par = GenMod_Par()  

  for seed in 1:n_seeds
    Random.seed!(seed)
    print(seed)

    for (i, perc_increase) in enumerate(percent_increases)
        #Modify attack rates by the percentage increase
        par = GenMod_Par()
        
        par.a_CR_1 = (par.a_CR_1 + rmax_high)
        par.a_CR_2 = (par.a_CR_2 + rmax_high)
        par.a_PC_1 = (par.a_PC_1 + rmax_high)
        par.a_PC_2 = (par.a_PC_2 + rmax_high)

        par.a_CR_1 *= (1 + perc_increase)
        par.a_CR_2 *= (1 + perc_increase)
        par.a_PC_1 *= (1 + perc_increase)
        par.a_PC_2 *= (1 + perc_increase)

        u0 = [1.0, 1.0, 0.5, 0.5, 0.25]
        prob_stoch = SDEProblem(GenMod_model!, GenMod_stochmodel!, u0, t_span, par)
        sol_stoch = solve(prob_stoch, reltol = 1e-15)
        grid_sol = sol_stoch(ts)

        #Record the mean, sd, and calculate CV
        maxhold[i, 1] = percent_increases[i]
        stdhold[i, 1] = std(grid_sol[5, 1500:2000])
        meanhold[i, 1] = mean(grid_sol[5, 1500:2000])
        cv_high_all[i, seed] = stdhold[i, 1] / meanhold[i, 1]

        #Reset parameters to run rmax low
        par = GenMod_Par()

        par.a_CR_1 = (par.a_CR_1 + rmax_low)
        par.a_CR_2 = (par.a_CR_2 + rmax_low)
        par.a_PC_1 = (par.a_PC_1 + rmax_low)
        par.a_PC_2 = (par.a_PC_2 + rmax_low)

        par.a_CR_1 *= (1 + perc_increase)
        par.a_CR_2 *= (1 + perc_increase)
        par.a_PC_1 *= (1 + perc_increase)
        par.a_PC_2 *= (1 + perc_increase)

        u0 = [1.0, 1.0, 0.5, 0.5, 0.25]
        prob_stoch = SDEProblem(GenMod_model!, GenMod_stochmodel!, u0, t_span, par)
        sol_stoch = solve(prob_stoch, reltol = 1e-15)
        grid_sol = sol_stoch(ts)

        # Recording the statistics
        maxhold[i, 1] = percent_increases[i]
        stdhold[i, 2] = std(grid_sol[5, 1800:2000])
        meanhold[i, 2] = mean(grid_sol[5, 1800:2000])
        cv_low_all[i, seed] = stdhold[i, 2] / meanhold[i, 2]
        
    end
end  

    ##Calculate change in percent CV from initial condition (aka 0)
    cv_high_avg = mean(cv_high_all, dims=2)
    cv_low_avg = mean(cv_low_all, dims=2)
    
    cv_rmax_high = cv_high_avg[1]
    percent_change_cv_high = [(cv - cv_rmax_high) / cv_rmax_high * 100 for cv in cv_high_avg]

    cv_rmax_low = cv_low_avg[1]
    percent_change_cv_low = [(cv - cv_rmax_low) / cv_rmax_low * 100 for cv in cv_low_avg]

    #Create a DataFrame
    df = DataFrame(
        Percent_Increase = vec(maxhold[:, 1]),
        Percent_Change_CV_High = vec(percent_change_cv_high[:, 1]), 
        Percent_Change_CV_Low = vec(percent_change_cv_low[:, 1]), 
    )

    #Write to CSV file for plotting in R
    #CSV.write("/Users/reillyoconnor/Desktop/Julia Projects/Slow-Fast/Julia/Theory GC Percent.csv", df)

    return
end


/