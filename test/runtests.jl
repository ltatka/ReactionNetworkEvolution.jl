using NetEvolve
using Test

@testset "Basic functions" begin

	r1 = NetEvolve.Reaction(["S2", "S1"], ["S3", "S4"], 5.0)
	r2 = NetEvolve.Reaction(["S1", "S2"], ["S4", "S3"], 1.3)
	myset = Set([r1, r2])
	@test length(myset) == 1


	settings, objfunct = NetEvolve.read_usersettings("DEFAULT")

	network1_str = """// Created by libAntimony v2.12.0
	// Compartments and Species:
	species S0, S1, S2;
	
	// Reactions:
	_J0: S0 -> S0 + S0; k1*S0;
	_J1: S1 -> S0; k2*S1;
	_J2: S1 + S2 -> S1 + S2; k3*S1*S2;
	_J3: S1 + S1 -> S0; k4*S1*S1;
	_J4: S1 + S2 -> S1 + S1; k5*S1*S2;
	_J5: S0 -> S0; k6*S0;
	_J6: S2 -> S1 + S2; k7*S2;
	_J7: S2 -> S1 + S1; k8*S2;
	_J8: S2 -> S2; k9*S2;
	_J9: S0 + S1 -> S2; k10*S0*S1;
	
	// Species initializations:
	S0 = 1;
	S1 = 5;
	S2 = 9;
	
	// Variable initializations:
	k1 = 33.3878793205463;
	k2 = 46.5172911588653;
	k3 = 52.5686815341686;
	k4 = 8.04051720757499;
	k5 = 51.9900520508139;
	k6 = 65.3963811855709;
	k7 = 0.983974104770924;
	k8 = 21.6057091196114;
	k9 = 7.05280598345276;
	k10 = 34.6814381533787;
	"""

	network2_str = """// Created by libAntimony v2.12.0
	// Compartments and Species:
	species S0, S1, S2;
	
	// Reactions:
	_J9: S0 + S1 -> S2; k1*S0*S1;
	_J0: S0 -> S0 + S0; k2*S0;
	_J1: S1 -> S0; k3*S1;
	_J3: S1 + S1 -> S0; k4*S1*S1;
	_J5: S0 -> S0; k5*S0;
	_J6: S2 -> S1 + S2; k6*S2;
	_J7: S2 -> S1 + S1; k7*S2;
	_J4: S1 + S2 -> S1 + S1; k8*S1*S2;
	_J8: S2 -> S2; k9*S2;
	_J2: S1 + S2 -> S1 + S2; k10*S1*S2;
	
	
	// Species initializations:
	S0 = 1;
	S1 = 5;
	S2 = 9;
	
	// Variable initializations:
	k1 = 33.3878793205463;
	k2 = 46.5172911588653;
	k3 = 52.5686815341686;
	k4 = 8.04051720757499;
	k5 = 51.9900520508139;
	k6 = 65.3963811855709;
	k7 = 0.983974104770924;
	k8 = 21.6057091196114;
	k9 = 7.05280598345276;
	k10 = 34.6814381533787;

	"""
	network1 = NetEvolve.convert_from_antimony(network1_str)
	@test [["S0"], ["S0", "S0"]] ∈ keys(network1.reactionlist)
	@test network1.reactionlist[[["S0"], ["S0", "S0"]]].rateconstant == 33.3878793205463

	network2 = NetEvolve.convert_from_antimony(network2_str)
	myset2 = Set([network1, network2])
	@test length(myset2) == 1

	ng = NetEvolve.get_networkgenerator(settings)
	@test ng.numreactions == 5

	population = NetEvolve.generate_network_population(settings, ng)
	println(length(population) == 100)
	@test typeof(population) == Vector{NetEvolve.ReactionNetwork}

	
end

@testset "solver tests" begin

	settings_dict = Dict{String, Any}(
    "specieslist" => ["A", "B"],
    "initialconditions" => [10, 0],
    )
	settings = NetEvolve.read_usersettings(settings_dict)
	smallnetwork_str = """A -> B; k1*A 
	k1 = 2
	A = 10
	B = 0"""
	smallnetwork = NetEvolve.convert_from_antimony(smallnetwork_str)
	tspan = (0.0, 5.0)
	sol = NetEvolve.solve_ode(smallnetwork, tspan)
	# Get rate of change for A at t = 0
	
	rateofchange = (sol.u[2][1] - sol.u[1][1])/sol.t[2]
	println(rateofchange)
	@test abs(rateofchange + 20) < 0.00001

	settings, objfunct = NetEvolve.read_usersettings("DEFAULT")
	astr="""
	S0 + S1 -> S2; k1*S0*S1
	S2 + S2 -> S2; k2*S2*S2
	S1 -> S0; k3*S1
	S2 -> S1 + S0; k4*S2

	k1 = 1
	k2 = 3
	k3 = 0.5
	k4 = 1

	S0 = 1
	S1 = 5
	S2 = 9
	"""
	network = NetEvolve.convert_from_antimony(astr)
	# Check initial rates of change
	du = NetEvolve.test_ode_funct([0.,0.,0.], [1., 5., 9.], network, 0.1)
	println(du)
	@test du == [6.5, 1.5, -13]
end

@testset "Network mutations" begin
	settings, objfunct = NetEvolve.read_usersettings("DEFAULT")
	astr="""
	S0 + S1 -> S2; k1*S0*S1
	S2 + S2 -> S2; k2*S2*S2
	S1 -> S0; k3*S1
	S2 -> S1 + S0; k4*S2

	k1 = 1
	k2 = 3
	k3 = 0.5
	k4 = 1

	S0 = 1
	S1 = 5
	S2 = 9
	"""
	network = NetEvolve.convert_from_antimony(astr)
	network = NetEvolve.mutaterateconstant(network, settings)
	keyslist = [[["S0", "S1"], ["S2"]],
				[["S2", "S2"], ["S2"]],
				[["S1"], ["S0"]],
				[["S2"], ["S0", "S1"]]]
	k_list = [1, 3, 0.5, 1]
	ismutated = false
	for (i, key) in enumerate(keyslist)
		if network.reactionlist[key].rateconstant != k_list[i]
			ismutated = true
		end
	end
	@test ismutated
	ng = NetEvolve.get_networkgenerator(settings)
	network = NetEvolve.addreaction(ng, network)
	@test length(keys(network.reactionlist)) == 5

	network = NetEvolve.deletereaction(network)
	@test length(keys(network.reactionlist)) == 5
	inactive_count = 0
	for key in keys(network.reactionlist)
		if !network.reactionlist[key].isactive
			inactive_count += 1
		end
	end
	@test inactive_count == 1
			
end

