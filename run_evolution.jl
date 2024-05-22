using ReactionNetworkEvolution
using ArgParse

s = ArgParseSettings()

@add_arg_table s begin
    "--ngenerations"
        arg_type = Int
        default = -1
    "--ntrials"
        arg_type = Int
        default = -1
    "--population_size"
        arg_type = Int
        default = -1
    "--pathtosettings"
        arg_type = String
        default = ""
    "--outputpath"
        arg_type = String
        default=""
    "--seed"
        arg_type = Int
        default = -1
    "--note"
        arg_type = String
        default = ""
end

parsed_args = parse_args(ARGS, s)

ReactionNetworkEvolution.run_evolution(ngenerations=parsed_args["ngenerations"],
              ntrials=parsed_args["ntrials"],
              population_size=parsed_args["population_size"],
              pathtosettings=parsed_args["pathtosettings"],
              outputpath=parsed_args["outputpath"],
              seed=parsed_args["seed"],
              note=parsed_args["note"])
