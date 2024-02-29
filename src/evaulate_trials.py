import os
import tellurium as te

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def load_model(path):
    f = open(path, "r")
    astr = f.read()
    f.close()
    r = te.loada(astr)
    return r

def writeout_model(r, path):
    f = open(path, "w")
    f.write(r.getCurrentAntimony())
    f.close()
def is_eigen_oscillator(r):
    try:
        r.steadyState()
        eigens = r.getFullEigenValues()
        for num in eigens:
            if num.real > 0 and num.imag != 0:
                return True
        return False
    except Exception as e:
        return False

def evaluate_models_dir(path, outputpath, bestonly=False):
    # This is for when you just have all the models in a directory
    models_evaluated = 0
    oscillators = 0
    for model_file in os.listdir(path):
        if bestonly:
            if model_file.startswith("bestmodel"):
                models_evaluated += 1
                full_model_path = os.path.join(path, model_file)
                r = load_model(full_model_path)
                if is_eigen_oscillator(r):
                    oscillators += 1
                    outpath = os.path.join(outputpath, model_file)
                    r.reset()
                    f = open(outpath, "w")
                    f.write(r.getCurrentAntimony())
                    f.close()
        else:
            models_evaluated += 1
            full_model_path = os.path.join(path, model_file)
            r = load_model(full_model_path)
            if is_eigen_oscillator(r):
                oscillators += 1
                outpath = os.path.join(outputpath, model_file)
                if os.path.isfile(outpath):
                    print("duplicate model ID found")
                    outpath += "_"
                f = open(outpath, "w")
                f.write(r.getCurrentAntimony())
                f.close()
    return models_evaluated, oscillators

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
        total, oscs = evaluate_models_dir(fullpath, outputpath, bestonly=bestonly)
        models_evaluated += total
        oscillators += oscs
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


inputpath = "/home/hellsbells/Desktop/networkEv/Data/BenchmarkTests/1000_trials"
outputpath = "/home/hellsbells/Desktop/networkEv/Data/BenchmarkTests/1000_trials_oscillators"

plot_trials(outputpath, 1,11)

# evaluate_trials_best_models(inputpath, outputpath, bestonly=True)

# inpath = "/home/hellsbells/Desktop/networkEv/Data/oscillators_from_lab_update"
# outpath = "/home/hellsbells/Desktop/networkEv/Data/oscillators_from_lab_update_true"
#
# # evaluate_trials_best_models(inpath, outpath, childdir=None)
#
# evaluate_models_dir(inpath, outpath, bestonly=False)