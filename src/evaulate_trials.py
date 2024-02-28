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

def is_eigen_oscillator(r):
    try:
        r.simulate(0, 1000, 5000)
        eigens = r.getFullEigenValues()
        for num in eigens:
            if num.real > 0 and num.imag != 0:
                return True
        return False
    except Exception as e:
        return False

def evaluate_trials_best_models(inputpath, outputpath, limit=None):
    if not os.path.isdir(outputpath):
        os.makedirs(outputpath)
    models_evaluated = 0
    oscillators = 0
    for file in os.listdir(inputpath):
        fullpath = os.path.join(inputpath, file)
        fullpath = os.path.join(fullpath, "final_models")
        for model_file in os.listdir(fullpath):
            if model_file.startswith("bestmodel"):
                models_evaluated += 1
                full_model_path = os.path.join(fullpath, model_file)
                r = load_model(full_model_path)
                if is_eigen_oscillator(r):
                    oscillators += 1
                    outpath = os.path.join(outputpath, model_file)
                    f = open(outpath, "w")
                    f.write(r.getCurrentAntimony())
                    f.close()
                if limit and models_evaluated >= limit:
                    print(f"Evaluated {models_evaluated} models and found {oscillators} oscillators")
                    print(f"Success rate: {oscillators / models_evaluated}")
                    return None
    print(f"Evaluated {models_evaluated} models and found {oscillators} oscillators")
    print(f"Success rate: {oscillators/models_evaluated}")


inputpath = "/home/hellsbells/Desktop/networkEv/Data/BenchmarkTests/1000_trials"
outputpath = "/home/hellsbells/Desktop/networkEv/Data/BenchmarkTests/1000_trials_oscillators"

evaluate_trials_best_models(inputpath, outputpath)