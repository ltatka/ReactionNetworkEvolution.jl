import os
import tellurium as te
from pathlib import Path
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

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
        if num.real > 0 and num.imag != 0:
            return True
    return False

def is_eigen_oscillator(r, writeoutpath):
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
                f = open(f"{writeoutpath}_flagged", "w")
                astr = r.getCurrentAntimony()
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
                if is_eigen_oscillator(r, outpath):
                    oscillators += 1
                    writeout_model(r, fitness, outpath)
                break
        else:
            full_model_path = os.path.join(path, model_file)
            r, fitness = load_model(full_model_path)
            models_evaluated += 1
            outpath = os.path.join(outputpath, model_file)
            if is_eigen_oscillator(r, outpath):
                oscillators += 1
                writeout_model(r, fitness, outpath)
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

astr = """
S0 -> S0 + S0; k1*S0
S0 + S1 -> S1; k2*S0*S1
S2 + S2 -> S1 + S2; k3*S2*S2
S0 + S2 -> S1 + S1; k4*S0*S2
S0 + S1 -> S1 + S1; k5*S0*S1
S2 -> S1 + S1; k6*S2
S2 -> S2 + S2; k7*S2
S0 + S2 -> S1 + S2; k8*S0*S2
S0 + S2 -> S0 + S2; k9*S0*S2
S0 -> S0 + S1; k10*S0
S2 -> S2; k11*S2
S0 + S2 -> S2; k12*S0*S2
S1 -> S0 + S0; k13*S1
S0 + S1 -> S0 + S0; k14*S0*S1
S2 -> S1; k15*S2
S1 -> S1; k16*S1
k1 = 2.272385348370609
k2 = 1.9020192182903577
k3 = 5.657195202706555
k4 = 16.639776452983458
k5 = 2.429586518569796
k6 = 1.8258797353770664
k7 = 34.060224078166655
k8 = 12.24006942737817
k9 = 2.4485979532961055
k10 = 18.86476105832696
k11 = 4.553327991133054
k12 = 14.912567209288481
k13 = 2.87090423417964
k14 = 6.142179542580343
k15 = 0.3809915986054561
k16 = 2.6863886555898477
S0 = 1.0
S1 = 5.0
S2 = 9.0

#fitness: 0.009997643507665454"""


inputpath = "/home/hellsbells/Desktop/networkEv/Data/BenchmarkTests/1000_trials_fixed"
outputpath = os.path.join(Path(inputpath).parent.absolute(), "1000_trials_fixed_oscillators_newalgo_fixed2")

evaluate_trials_best_models(inputpath, outputpath, bestonly=True)

