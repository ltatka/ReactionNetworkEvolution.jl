using RoadRunner
using Random

include("network_cleanup.jl")
include("reaction_network.jl")

function has_oscillator_eigens(eigen_array)
    for i in 1:size(eigen_array)[1]
        real = eigen_array[i, 1]
        imag = eigen_array[i, 2]
        if real >= 0 && imag != 0
            return true
        end
    end
    return false
end

function is_oscillator(astr::String)
    r = RoadRunner.loada(astr)
    try
        RoadRunner.steadyState(r) # WTF
        eigens = RoadRunner.getEigenvalues(r)
        concentrations = RoadRunner.getFloatingSpeciesConcentrations(r)
        
        # If it has the correct eigens and positive concentrations, return true
        return has_oscillator_eigens(eigens) && all(>=(0), concentrations)
    catch 
        # Sometimes it will error when calculating the eigens or steady state, but that doesn't necessarily mean
        # it is not an oscillator
    end
    # If it fails to get eigens values, try to simulate it
    # If the simulation fails, adjust tolerances
    try
        RoadRunner.resetToOrigin(r)
        RoadRunner.setTimeStart(r, 0)
        RoadRunner.setTimeEnd(r, 50)
        RoadRunner.setNumPoints(r, 100000)
        RoadRunner.simulate(r)
    catch
        try # Try to simulate with adjusted tolerance
            RoadRunner.resetToOrigin(r)
            RoadRunner.setTimeStart(r, 0)
            RoadRunner.setTimeEnd(r, 50)
            RoadRunner.setNumPoints(r, 100000)
            RoadRunner.setCurrentIntegratorParameterString(r, "relative_tolerance", 1e-10)
            RoadRunner.simulate(r)
        catch # if that doesn't work either, it's not an oscillator
            return false
        end
    end
    # Try to get the eigens again, if that doesn't work, it's not an oscillator
    try
        eigens = RoadRunner.getEigenvalues(r)
        #Now check if the eigen values indicate an oscillator
        # if it does have oscillator eigens, it's an oscillator
        if has_oscillator_eigens(eigens)
            return true
        end
    catch 
        return false
    end
    # 
    # If that still didn't work, calculate steady state again
    try
        RoadRunner.steadyState(r)
        eigens = RoadRunner.getEigenvalues(r)
        concentrations = RoadRunner.getFloatingSpeciesConcentrations(r)
        return has_oscillator_eigens(eigens) && all(>=(0), concentrations)
    catch
        # If after all that, something still doesn't work, then it's not an oscillator
        return false
    end
end

function dummy_is_oscillator(astr::String)
    return rand() < 0.5
end

function make_output_dirs(outputpath::String)
    if !isdir(joinpath(outputpath, "SUCCESS"))
        mkdir(joinpath(outputpath, "SUCCESS"))
    end
    if !isdir(joinpath(outputpath, "FAIL"))
        mkdir(joinpath(outputpath, "FAIL"))
    end
end

function process_oscillators(outputpath::String)
    make_output_dirs(outputpath)
    for directory in readdir(outputpath)
        if directory != "SUCCESS" && directory != "FAIL" # Ignore output folders
            # Navigate to the inner directories that contain final_models, datatracker.json, and settings.json
            for model_directory in readdir(joinpath(outputpath, directory))
                # Check if there's a final_models directory
                if !isdir(joinpath(outputpath, directory, model_directory, "final_models"))
                    continue
                end
                # If the directory does exist, get the best model from model_directory
                for model_file in readdir(joinpath(outputpath, directory, model_directory, "final_models"))
                    if !startswith(model_file, "bestmodel")
                        continue
                    end
                    astr = load_antimony_file(joinpath(outputpath, directory, model_directory, "final_models", model_file))
                    is_oscillator_model = dummy_is_oscillator(astr)
                    if is_oscillator_model
                        destination = joinpath(outputpath, "SUCCESS", model_directory)
                    else
                        destination = joinpath(outputpath, "FAIL", model_directory)
                    end
                    mv(joinpath(outputpath, directory, model_directory), destination)
                    break
                end
            end
        # Delete the (now) empty directory
        rm(joinpath(outputpath, directory))
        end
    end
end

outputpath = "/home/hellsbells/Desktop/evolution_output/test_sort"

process_oscillators(outputpath)





# astr = """// Created by libAntimony v2.12.0
# // Compartments and Species:
# species S0, S1, S2;

# // Reactions:
# _J0: S0 -> S0 + S0; k1*S0;
# _J1: S0 + S1 -> S1; k2*S0*S1;
# _J2: S0 + S1 -> S1 + S1; k3*S0*S1;
# _J3: S0 + S2 -> S0 + S0; k4*S0*S2;
# _J4: S1 -> S0 + S1; k5*S1;
# _J5: S0 -> S2 + S2; k6*S0;
# _J6: S2 -> S2; k7*S2;
# _J7: S1 + S2 -> S2; k8*S1*S2;
# _J8: S2 -> S1; k9*S2;

# // Species initializations:
# S0 = 1;
# S1 = 5;
# S2 = 9;

# // Variable initializations:
# k1 = 82.0537175603692;
# k2 = 1.46247635550851;
# k3 = 11.0307336282238;
# k4 = 71.8330483162337;
# k5 = 0.500520887516744;
# k6 = 82.8386889676032;
# k7 = 31.2822563337994;
# k8 = 25.2435734826658;
# k9 = 9.24024972949223;

# // Other declarations:
# const k1, k2, k3, k4, k5, k6, k7, k8, k9;
# """

# println(is_oscillator(astr))