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