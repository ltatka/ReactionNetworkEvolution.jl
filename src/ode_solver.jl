using DifferentialEquations
using Sundials
# using Plots
# using BenchmarkTools

include("reaction_network.jl")




function ode_funct!(du::Array{Float64}, u::Array{Float64}, network::ReactionNetwork, t::Float64)
    specieslist = network.specieslist
    
    # Reset du
    for i = 1:length(specieslist)
        du[i] = 0.0
    end

    for key in keys(network.reactionlist)
        reaction = network.reactionlist[key]
        if !reaction.isactive
            continue
        end
        # Get the relevant concentrations
        dspecies = 1 # Is there a case where this would be wrong?
        for s in reaction.substrate
            idx = findfirst(item -> item == s, specieslist)
            dspecies *= u[idx]
        end
        # Multiply by the rate constant to get the rate for *this* reaction
        dspecies *= reaction.rateconstant
        # Subtract this rate for substrates
        for s in reaction.substrate
            idx = findfirst(item -> item == s, specieslist)
            du[idx] -= dspecies
        end
        # Add this rate for products
        for p in reaction.product
            idx = findfirst(item -> item == p, specieslist)
            du[idx] += dspecies
        end
    end
    # for boundary species, reset the rate of change to 0
    for s in network.boundaryspecies
        idx = findfirst(item -> item == s, specieslist)
        du[idx] = 0.0
    end
end

function solve_ode(objfunct::ObjectiveFunction, network::ReactionNetwork)
    # 
    #Get time info
    t0 = first(objfunct.time)
    t_end = last(objfunct.time)
    stepsize = objfunct.time[2] - t0 # for now we assume even spacing of time points
    tspan = (t0, t_end)

    u0 = network.initialcondition
    ode_prob = ODEProblem(ode_funct!, u0, tspan, network)
    sol = solve(ode_prob, CVODE_BDF(), saveat=stepsize, verbose=false)
    return sol
    
end


function evaluate_fitness(objfunct:: ObjectiveFunction, network::ReactionNetwork)
    try 
        sol = solve_ode(objfunct, network)
        if length(sol.t) != length(objfunct.time) # If the time points are unequal, then simulation has failed.
            return DEFAULT_FITNESS
        end
        fitness = 0.0
        nspecies = length(network.floatingspecies)
        for i in 1:nspecies
            timeseries = sol[i,:]
            tempfitness = 0.0
            for (j, pt) in enumerate(timeseries)
                tempfitness += abs(pt - objfunct.objectivedata[j])
            end
            if 1/tempfitness > fitness
                fitness = 1/tempfitness
            end
        end
        return fitness
    catch e
        return DEFAULT_FITNESS
    end
    
end





