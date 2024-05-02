using NetEvolve
using ArgParse

s = ArgParseSettings()

@add_arg_table s begin
    "--ngenerations"
        arg_type = Int
        default = -1
    "--ntrials"
        arg_type = Int
        default = -1
    "--populationsize"
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

NetEvolve.run_evolution(ngenerations=parsed_args["ngenerations"],
              nbatches=parsed_args["ntrials"],
              populationsize=parsed_args["populationsize"],
              pathtosettings=parsed_args["pathtosettings"],
              outputpath=parsed_args["outputpath"],
              seed=parsed_args["seed"],
              note=parsed_args["note"])
