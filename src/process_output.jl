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

function fix_broken_oscillator(astr::String)
    # Sometimes models that have oscillator Eigen values but negative steady state concentrations
    # can be fixed by removing a single reaction
    astr_lines = split(astr, "\n") 
    for (i, line) in enumerate(astr_lines)
        # comment out reaction
        if occursin("->", line)
            astr_lines[i] = "#" * astr_lines[i]
            new_astr = join(astr_lines, "\n")
            if is_oscillator(new_astr)
                return new_astr #???
            end
            # If that didn't work, change it back
            astr_lines[i] = astr_lines[i][2:end]
        end
    end
    return "FAIL"
end

function is_broken_oscillator(astr::String)
    # Sometimes models that have oscillator Eigen values but negative steady state concentrations
    # can be fixed by removing a single reaction
    r = RoadRunner.loada(astr)
    try
        RoadRunner.steadyState(r) # WTF
        eigens = RoadRunner.getEigenvalues(r)
        concentrations = RoadRunner.getFloatingSpeciesConcentrations(r)
        # If it has the correct eigens and positive concentrations, return true
        return has_oscillator_eigens(eigens) && !all(>=(0), concentrations)
    catch 
    end
    return false
end

function is_oscillator(astr::String)
    r = RoadRunner.loada(astr)
    try
        RoadRunner.steadyState(r) 
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

function make_output_dirs(outputpath::String)
    timestamp = format(now(), "YYYYmmdd_HHMMSS")
    results_parent_dir = joinpath(outputpath, "results_$timestamp")
    if !isdir(results_parent_dir)
        mkdir(results_parent_dir)
    end
    if !isdir(joinpath(results_parent_dir, "SUCCESS"))
        mkdir(joinpath(results_parent_dir, "SUCCESS"))
    end
    if !isdir(joinpath(results_parent_dir, "FAIL"))
        mkdir(joinpath(results_parent_dir, "FAIL"))
    end
    return results_parent_dir
end

function process_oscillators(outputpath::String)
    results_parent_dir = make_output_dirs(outputpath)
    
    for directory in readdir(outputpath)
        if !startswith(directory, "results") # Ignore results folders
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
                    is_oscillator_model = is_oscillator(astr)
                    # If it's an oscillator, move it to the SUCCESS dir
                    if is_oscillator_model
                        destination = joinpath(results_parent_dir, "SUCCESS", model_directory)
                    else
                        # Try to fix it 
                        if is_broken_oscillator(astr)
                            astr = fix_broken_oscillator(astr)
                            if astr != "FAIL" # If fixing it worked, overwrite file
                                model_file_path = joinpath(outputpath, directory, model_directory, "final_models", model_file)
                                rm(model_file_path) # Remove the old filepath
                                open(model_file_path, "w") do file
                                    write(file, astr)    
                                end
                                # Fixed file will be moved to the successful results directory
                                destination = joinpath(results_parent_dir, "SUCCESS", model_directory)
                            end
                        end
                        destination = joinpath(results_parent_dir, "FAIL", model_directory)
                    end
                    mv(joinpath(outputpath, directory, model_directory), destination)
                    break
                end
            end
            
            # Delete the (now) empty directory
            try
                rm(joinpath(outputpath, directory))
            catch e
                println(e)
            end
        end
    end
    println("Output written to $results_parent_dir")
end