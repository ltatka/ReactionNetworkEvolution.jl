import JSON
include("reaction_network.jl")
include("evo_utils.jl")
include("settings.jl")
include("ode_solver.jl")
include("network_cleanup.jl")
using Plots


function plot_timeseries_from_file(objfunct:: ObjectiveFunction, network_path::String, save_path::String)
    network = convert_from_antimony(network_path)
    solution = solve_ode(objfunct, network)
    plt = Plots.plot(solution)
    # save_pa = dirname(pwd()) * "$(network.ID).png"
    savefig(plt, save_path)
    return plt
end


inpath = "/home/hellsbells/Desktop/networkEv/Data/BenchmarkTests/1000_trials_fixed_oscillators_all/bestmodel_n7QEsETEzFmy"
outpath = "/home/hellsbells/Desktop/networkEv/Data/BenchmarkTests/1000_trials_fixed_oscillators_all/test"
pathtosettings = "/home/hellsbells/Desktop/networkEv/test_files/updownObjFunc.json"

settings = read_usersettings(pathtosettings)
objfunct = get_objectivefunction(settings)

network = convert_from_antimony(inpath)
f = evaluate_fitness(objfunct, network)
println(f)
# plot_timeseries_from_file(objfunct, inpath, outpath)

# rs = ReactionSet()
# pushreaction(rs, r1)
# pushreaction(rs, r2)
# display(rs)

# rs2 = ReactionSet()
# pushreaction(rs2, r1, addrateconstants=false)
# pushreaction(rs2, r2, addrateconstants=false)
# display(rs2)

# mkdir("Thisisatest")


# model = """ # This is the model that timeseries test file comes from
# S2 -> S2 + S2; k1*S2
# S2 + S3 -> S1; k2*S2*S3
# S2 -> S1 + S3; k3*S2
# S1 + S1 -> S1; k4*S1*S1
# S2 -> S1; k5*S2
# k1 = 0.3945790331714698
# k2 = 0.9318069012395747
# k3 = 0.9759098143496107
# k4 = 0.46279808093283303
# k5 = 1.2774098437413
# S1 = 1.0
# S2 = 2.0
# S3 = 2.0"""



# filepath = "/home/hellsbells/Desktop/test_ant.txt"

# network = convert_from_antimony(filepath)

# display(network)

# specieslist = network.specieslist
# initialcondition = network.initial_conditions
# numreactions = 6


# ng = NetworkGenerator(specieslist, initialcondition, numreactions, reactionprobabilities, rateconstantrange)


macro timeout(seconds, expr, fail)
    quote
        tsk = @task $expr
        schedule(tsk)
        Timer($seconds) do timer
            istaskdone(tsk) || Base.throwto(tsk, InterruptException())
        end
        try
            fetch(tsk)
        catch _
            $fail
        end
    end
end

x = @timeout 1 begin
    sleep(0.1)
    println("done")
    1
end "failed"

println(x)

println("Done")