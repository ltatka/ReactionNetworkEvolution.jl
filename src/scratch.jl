import JSON
include("reaction_network.jl")

# s = """{
# "test": 2}
# """

# j = JSON.parse(s)

# println(j["test"])


# j = JSON.parsefile("/home/hellsbells/Desktop/networkEv/src/test.json")
# println(j["species"])

# println(typeof(j))
# println(keys(j))

# a = Dict("thing" => "value", "otherthing" => 3)

# struct Test
#     thing1::Float64
#     thing2::Float64

#     function Test(p)
#         new(p[1], p[2])
#     end

#     function Test(p1, p2)
#         new(p1, p2)
#     end

# end

# t = Test([1,2])

# display(t)
# t = Test(1, 2)
# display(t)

r1 = Reaction(["S1", "S2"], ["S3"], .25)
r2 = Reaction(["S2", "S1"], ["S3"], .35)

print(reaction_isequal(r1, r2))

struct ReactionSet
    reactionlist:: Vector{Reaction}
    
    function ReactionSet()
        new([])
    end

end

function pushreaction(reactionset::ReactionSet, reaction; addrateconstants=true)
    for r in reactionset.reactionlist
        if reaction_isequal(r, reaction)
            if addrateconstants
                r.rateconstant += reaction.rateconstant
            end
            return reactionset
        end
    end
    push!(reactionset.reactionlist, reaction)
    return reactionset
end

# rs = ReactionSet()
# pushreaction(rs, r1)
# pushreaction(rs, r2)
# display(rs)

# rs2 = ReactionSet()
# pushreaction(rs2, r1, addrateconstants=false)
# pushreaction(rs2, r2, addrateconstants=false)
# display(rs2)

mkdir("Thisisatest")


model = """ # This is the model that timeseries test file comes from
S2 -> S2 + S2; k1*S2
S2 + S3 -> S1; k2*S2*S3
S2 -> S1 + S3; k3*S2
S1 + S1 -> S1; k4*S1*S1
S2 -> S1; k5*S2
k1 = 0.3945790331714698
k2 = 0.9318069012395747
k3 = 0.9759098143496107
k4 = 0.46279808093283303
k5 = 1.2774098437413
S1 = 1.0
S2 = 2.0
S3 = 2.0"""