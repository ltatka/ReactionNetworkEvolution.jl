import JSON

s = """{
"test": 2}
"""

j = JSON.parse(s)

println(j["test"])


j = JSON.parsefile("/home/hellsbells/Desktop/networkEv/src/test.json")
println(j["species"])

println(typeof(j))
println(keys(j))

a = Dict("thing" => "value", "otherthing" => 3)

struct Test
    thing1::Float64
    thing2::Float64

    function Test(p)
        new(p[1], p[2])
    end

    function Test(p1, p2)
        new(p1, p2)
    end

end

t = Test([1,2])

display(t)
t = Test(1, 2)
display(t)