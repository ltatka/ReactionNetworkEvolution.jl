using networkEv
using Test

@testset "networkEv.jl" begin
	 
	r1 = networkEv.Reaction(["S2", "S1"], ["S3", "S4"], 5.0)
	r2 = networkEv.Reaction(["S1", "S2"], ["S4", "S3"], 1.3)
	myset = Set([r1, r2])
	@test length(myset) == 1

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
	
	// Other declarations:
	const k1, k2, k3, k4, k5, k6, k7, k8, k9, k10;"""

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
	
	// Other declarations:
	const k1, k2, k3, k4, k5, k6, k7, k8, k9, k10;
	"""
	network1 = convert_from_antimony_string(network1_str)
	network2 = convert_from_antimony_string(network2_str)
	myset2 = Set([network1, network2])
	@test length(myset2) == 1

end
