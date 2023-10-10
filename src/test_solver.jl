using Plots
using CSV
using DataFrames
using JSON
include("ode_solver.jl")
include("evo_utils.jl")
# include("settings.jl")

pathtosettings = "/home/hellsbells/Desktop/networkEv/src/test.json"

settings = read_usersettings(pathtosettings)

println(settings.reactionprobabilities)


# # #### TEST
# usersettings = UserSettings(["S1", "S2", "S3"], 
#     [1.0, 2.0, 2.0], 
#     "/home/hellsbells/Desktop/testtimeseries.csv", 
#     ["S1","S2", "S3"],
#     10,
#     200)

objfunct = get_objectivefunction(settings)
ng = get_networkgenerator(settings)

# settings = Settings(.1, ReactionProbabilities(0.25, 0.25, 0.25, 0.25), MutationProbabilities(0.5, 0.5), usersettings)
# rn = NetworkGenerator(["S1", "S2", "S3"], [1.0, 2.0, 2.0], 5, ReactionProbabilities(.25, 0.25, 0.25, 0.25), [0.1, 2.0])
# # network = generate_random_network(rn)

# population = evolve(settings, rn, objfunct)

# bestnetwork = population[1]
# println(convert_antimony(bestnetwork))

# df = DataFrame(CSV.File(settings.objectivedatapath))

# time = df[!, "time"]
# objective = df[!, ["S2", "S1"]]
# print(typeof(objective))



# tspan = (0.0, 1)

# u0 = network.initialcondition
# ode_prob = ODEProblem(ode_funct!, u0, tspan, network)
# sol = solve(ode_prob, CVODE_BDF(), saveat=1.0/9)

# function solve_ode(network)
#     tspan = (0.0, 1)

#     u0 = network.initialcondition
#     ode_prob = ODEProblem(ode_funct!, u0, tspan, network)
#     sol = solve(ode_prob, CVODE_BDF(), saveat=1.0/9)
#     return sol
# end

# sol = solve_ode(objfunct, network)
# for (i,row) in enumerate(sol.u)
#     println(typeof(row[2]))
#     println(row[2])
# end

# fitness =evaluate_fitness(objfunct, network)
# print(fitness)
# # # plot(sol)
# # savefig("/home/hellsbells/Desktop/plot1.png")

# display(sol)
# display(sol.u)

# Get fitness
# fitness = 0.0
# for (j, row) in enumerate(sol.u)
#     global fitness
#     diff = abs(objective[j] - row[2])
#     fitness += diff
    
# end
# display(fitness)




println("Success")

