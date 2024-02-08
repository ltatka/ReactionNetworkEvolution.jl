using DifferentialEquations
using Sundials
using Plots
using BenchmarkTools
include("reaction_network.jl")

# # This is a macro that will kill a function if it takes too long. Useful 
# # for troublesome differential equations
# macro timeout(seconds, expr, fail)
#     quote
#         tsk = @task $expr
#         schedule(tsk)
#         Timer($seconds) do timer
#             istaskdone(tsk) || Base.throwto(tsk, InterruptException())
#         end
#         try
#             fetch(tsk)
#         catch _;
#             println("we got to the fail macro")
#             $fail
#         end
#     end
# end

struct ObjectiveFunction
    objectivespecies::Vector{String} # Can be more than 1
    objectivedata:: DataFrame 
    time::Vector{Float64}
    indexbyspecies::Dict
end

function get_objectivefunction(settings::Settings)

    objectivedata = DataFrame(CSV.File(settings.objectivedatapath))
    time = objectivedata[!, "time"]
    println("in get_objectivefunction: timepts: $(length(time))")
    indexbyspecies = Dict()
    for s in settings.objectivespecies
        idx = findfirst(item -> item == s, settings.specieslist)
        indexbyspecies[s] = idx
    end
    return ObjectiveFunction(settings.objectivespecies, objectivedata, time, indexbyspecies)
end

function ode_funct!(du, u, network::ReactionNetwork, t)
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

function solve_ode(objfunct, network)
    # 
    #Get time info
    t0 = first(objfunct.time)
    t_end = last(objfunct.time)
    stepsize = objfunct.time[2] - t0 # for now we assume even spacing of time points
    tspan = (t0, t_end)

    u0 = network.initialcondition
    ode_prob = ODEProblem(ode_funct!, u0, tspan, network)
    sol = solve(ode_prob, CVODE_BDF(), saveat=stepsize, verbose=false)
    # sol = @timeout MAX_TIME begin
    #     println(MAX_TIME)
    #     solve(ode_prob, CVODE_BDF(), saveat=stepsize, verbose=false)
    #     println("succss")
    # end println("fail")
    return sol
    
end


function evaluate_fitness(objfunct:: ObjectiveFunction, network::ReactionNetwork; sizepenalty=false)
    # This function evaluates the fitness of an individual based on how well the timeseries of the 
    # target species matches that of the target timeseries. 
    # Changing all fitness stuff here from now on
    TARGET_NETWORK_SIZE = 5 # Trying to encourage networks to be this size
    PENALTY_MULTIPLIER = 0.95
    try
        sol = solve_ode(objfunct, network)
        if length(sol.t) != length(objfunct.time) # If the time points are unequal, then simulation has failed.
            return DEFAULT_FITNESS
        end
        fitness = 0.0
        for (i, row) in enumerate(sol.u)
            for s in objfunct.objectivespecies
                idx = objfunct.indexbyspecies[s]
                fitness += abs(objfunct.objectivedata[!, s][i] - row[idx])
            end
        end
        if sizepenalty
            ## TRYING TO PENALIZE LARGE NETWORKS 
            # Take off some percentage of the total fitness according to how large the network is?
            numreactions = length(keys(network.reactionlist))
            excessreactions = numreactions - TARGET_NETWORK_SIZE 
            if excessreactions < 0 # We don't care if the network is too small
                excessreactions = 0
            end
            penalty = (1 - excessreactions/numreactions)*PENALTY_MULTIPLIER
        else
            penalty = 1
        end
        fitness = 1/fitness # Make larger fitness = better
        fitness = fitness*penalty

        return fitness # Or should this also assign the fitness to the network?
    catch e
        println(e)
        return DEFAULT_FITNESS
    end
end



function plot_timeseries(objfunct:: ObjectiveFunction, network:: ReactionNetwork; path=nothing)
    solution = solve_ode(objfunct, network)
    plt = plot(solution)
    if isnothing(path)
        path = dirname(pwd()) * "$(network.ID).png"
    end
    savefig(plt, path)
    return plt
end

