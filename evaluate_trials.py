import os
import tellurium as te
from pathlib import Path
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import shutil

def load_model(path):
    f = open(path, "r")
    astr = f.read()
    f.close()
    r = te.loada(astr)
    for line in astr.split("\n"):
        if line.startswith("#fitness"):
            f = line.split(" ")[1]
            f = float(f)
    return r, f

def writeout_model(r, fitness, path):
    f = open(path, "w")
    r.reset()
    f.write(r.getCurrentAntimony() + f"\n#fitness: {fitness}")
    f.close()

def check_eigens(eigen_array):
    for num in eigen_array:
        if num.real >= 0 and num.imag != 0:
            return True
    return False

def is_eigen_oscillator(r, fitness=None, writeoutpath=None):
    hasgoodvals = False
    # 1. Measure eigen values
    # eigens = r.getFullEigenValues()
    # hasgoodvals = check_eigens(eigens)
    # if hasgoodvals:
    #     return True
    # 2. If that didn't work, measure steady state
    try:
        r.steadyState()
        eigens = r.getFullEigenValues()
        hasgoodvals = check_eigens(eigens)
        if hasgoodvals:
            if r.S0 < 0 or r.S1 < 0 or r.S2 < 0:
                if writeoutpath is not None:
                    f = open(f"{writeoutpath}_flagged", "w")
                    r.reset()
                    astr = r.getCurrentAntimony()
                    if fitness is not None:
                        astr += f"\n#fitness: {fitness}"
                    f.write(astr)
                    f.close()
                return False
            else:
                return True
    except:
        pass
    # 3. If that didn't work, simulate for a bit
    r.reset()
    try:
        r.simulate(0, 50, 100000)

    except:
        try:
            r.integrator.relative_tolerance = 1e-10
            r.reset()
            r.simulate(0, 50, 100000)
        except:
            return False
    try:
        eigens = r.getFullEigenValues()
    except:
        return False
    hasgoodvals = check_eigens(eigens)
    if hasgoodvals:
        return True
    # 4. If that still didn't work, calculate steady state again
    try:
        r.steadyState()
        eigens = r.getFullEigenValues()
        hasgoodvals = check_eigens(eigens)
    except:
        pass
    # And if that doesn't work then it's definitely (?) not an oscillator
    return hasgoodvals

def is_eigen_oscillator_old(r, f, writeoutpath):
    try:
        r.steadyState()
        eigens = r.getFullEigenValues()
        for num in eigens:
            if num.real > 0 and num.imag != 0:
                return True
        return False
    except Exception as e:
        try:
            r.simulate(0, 10000, 1000)
            r.steadyState()
            eigens = r.getFullEigenValues()
            for num in eigens:
                if num.real > 0 and num.imag != 0:
                    return True
            return False
        except Exception as e:
            writeout_model(r, f, writeoutpath)
            return False

def evaluate_models_dir(path, outputpath, bestonly=False):
    # This is for when you just have all the models in a directory
    models_evaluated = 0
    oscillators = 0
    for model_file in os.listdir(path):
        if bestonly:
            if model_file.startswith("bestmodel"):
                # print(path, model_file)
                models_evaluated += 1
                full_model_path = os.path.join(path, model_file)
                r, fitness = load_model(full_model_path)
                outpath = os.path.join(outputpath, model_file)
                if is_eigen_oscillator(r, fitness=fitness, writeoutpath=outpath):
                    oscillators += 1
                    writeout_model(r, fitness, outpath)
                break
        else:
            full_model_path = os.path.join(path, model_file)
            r, fitness = load_model(full_model_path)
            models_evaluated += 1
            outpath = os.path.join(outputpath, model_file)
            if is_eigen_oscillator(r, fitness=fitness, writeoutpath=outpath):
                oscillators += 1
                writeout_model(r, fitness, outpath)
    return models_evaluated, oscillators

def copy_(src, dst):
    if os.path.isdir(src):
        shutil.copytree(src, dst)
    else:
        shutil.copy(src, dst)
def evaluate_trials_best_models(inputpath, outputpath, childdir="final_models", limit=None, bestonly=False):
    print("starting...")
    if not os.path.isdir(outputpath):
        os.makedirs(outputpath)
    models_evaluated = 0
    oscillators = 0
    for file in os.listdir(inputpath): # Eg. 2024-02-24T12:03:00.332
        fullpath = os.path.join(inputpath, file)
        if childdir:
            fullpath = os.path.join(fullpath, childdir)
        try:
            total, oscs = evaluate_models_dir(fullpath, outputpath, bestonly=bestonly)
            models_evaluated += total
            oscillators += oscs
            if oscs > 0:
                copy_(inputpath, outputpath)
        except:
            continue
        if models_evaluated % 50 == 0:
            print(f"Evaluated: {models_evaluated}")
        # for model_file in os.listdir(fullpath):
        #     total, oscs = evaluate_models_dir(fullpath, outputpath, bestonly=bestonly)
        #     models_evaluated += total
        #     if models_evaluated % 50 == 0:
        #         print(f"Evaluated {models_evaluated}")
        #     oscillators += oscs

            # if limit and models_evaluated >= limit:
            #         print(f"Evaluated {models_evaluated} models and found {oscillators} oscillators")
            #         print(f"Success rate: {oscillators / models_evaluated}")
            #         return None
    print(f"Evaluated {models_evaluated} models and found {oscillators} oscillators")
    print(f"Success rate: {oscillators/models_evaluated}")
    return oscillators, models_evaluated

def plot_trials(inputpath, rows, cols, childdir="final_models"):
    limit = rows*cols
    count = 0
    fig, axs = plt.subplots(rows, cols)
    models = []
    for file in os.listdir(inputpath): # Eg. 2024-02-24T12:03:00.332
        fullpath = os.path.join(inputpath, file)
        if childdir:
            fullpath = os.path.join(fullpath, childdir)
        for model_file in os.listdir(fullpath):
            if model_file.startswith("bestmodel"):
                count += 1
                models.append(load_model(os.path.join(fullpath, model_file)))
                break
        if count > limit:
            break
    idx = 0
    for i in range(rows):
        for j in range(cols):
            try:
                r = models[idx]
                m = r.simulate()
                axs[i,j].plot(m['time'], m[:,1:])
                # a
            except:
                pass
            idx += 1
            if idx >= len(models):
                break
    plt.show()

def fix_model(astr):
    split_astr = astr.split("\n")
    for i in range(len(split_astr)):
        if "->" in split_astr[i]:
            split_astr[i] = "#" + split_astr[i]
            newastr = "\n".join(split_astr)
            r = te.loada(newastr)
            if is_eigen_oscillator(r):
                return newastr
            # If the model is not fixed, then uncomment the line
            split_astr[i] = split_astr[i][1:]
    # If the model can't be fixed, return None
    return None
def fix_flagged_models(inputpath):
    # input path is where all the final oscillators are stored.
    num_fixed = 0
    for file in os.listdir(inputpath):
        if "flagged" in file:
            full_model_path = os.path.join(inputpath, file)
            f = open(full_model_path, "r")
            astr = f.read()
            f.close()
            fitness = -1
            for line in astr.split("\n"):
                if line.startswith("#fitness"):
                    fitness = line.split(" ")[1]
                    fitness = float(fitness)
            fixed_astr = fix_model(astr)
            if fixed_astr is not None:
                fixed_astr += f"\n#fitness(old): {fitness}"
                f = open(f"{full_model_path}_fixed", "w")
                f.write(fixed_astr)
                f.close()
                print(f"successfully fixed model {file}")
                num_fixed += 1
            else:
                print("unable to fix model")
    return num_fixed



inputpath = "/home/hellsbells/Desktop/evolution_output/unique_species_no_spec/output"
outputpath = inputpath + "_success"

oscillators, total = evaluate_trials_best_models(inputpath, outputpath, bestonly=True)
print("attempting to fix flagged models")
fixed_models = fix_flagged_models(outputpath)
print(f"fixed {fixed_models} models\nnew success rate: {(oscillators+fixed_models)/total}")
print(f"Total: {total}\nFixed: {fixed_models}\nOriginal oscillators: {oscillators} ")

