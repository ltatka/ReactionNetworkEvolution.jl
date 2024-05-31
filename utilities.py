import tellurium as te
import math
import os
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
            fitness = line.split(" ")[1]
            fitness = float(fitness)
    return r, fitness

def check_eigens(eigen_array):
    for num in eigen_array:
        if num.real >= 0 and num.imag != 0:
            return True
    return False

def is_oscillator(r, fitness=None, writeoutpath=None):
    try:
        s = r.steadyState()
        eigens = r.getFullEigenValues()
        hasgoodvals = check_eigens(eigens)
        if hasgoodvals:
            if not all(i >= 0 for i in r.getFloatingSpeciesConcentrations()):
                if writeoutpath is not None:
                    f = open(f"{writeoutpath}_flagged", "w")
                    r.resetToOrigin()
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
    #r.reset()
    r.resetToOrigin()
    try:
        r.simulate(0, 50, 100000)

    except:
        try:
            r.integrator.relative_tolerance = 1e-10
            r.resetToOrigin()
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
        s = r.steadyState()
        change = r.getRatesOfChange()
        print(change)
        eigens = r.getFullEigenValues()
        hasgoodvals = check_eigens(eigens)
    except:
        pass
    # And if that doesn't work then it's definitely (?) not an oscillator
    return hasgoodvals and all(i >= 0 for i in r.getFloatingSpeciesConcentrations())


#--------------------------------------------------
# Plotting Timeseries Data
#-------------------------------------------------
def get_best_dimensions(n):
    # Find the best dimensions for a subplot grid to plot n subplots
    # Input: n, Int, the total number of items to plot
    # Returns: w, h, Ints, width and height for optimal arrangement
    sqrt = n**0.5
    closest_2 = [math.floor(sqrt)**2, math.ceil(sqrt)**2]
    differences = [abs(x - n) for x in closest_2]
    closest = closest_2[differences.index(min(differences))]
    best = int(closest**0.5)
    if closest >= n:
        return best, best
    else:
        return best, best + 1

def plot_timeseries_path(path, start, end, numpoints, savepath):
    plt.rcParams.update({'font.size': 6})
    filenames = os.listdir(path)
    n = len(filenames)
    rows, cols = get_best_dimensions(n)

    models = []
    for file in os.listdir(path):
        r, _ = load_model(os.path.join(path, file))
        models.append(r)
    idx = 0
    fig, axs = plt.subplots(rows, cols)
    for i in range(rows):
        for j in range(cols):
            r = models[idx]
            m = r.simulate(start, end, numpoints)
            axs[i, j].plot(m['time'], m[:, 1:])
            axs[i, j].set_title(filenames[idx])
            idx += 1
            if idx > n - 1:
                break
        if idx > n - 1:
            break
    fig.tight_layout()
    if savepath:
        plt.savefig(savepath)
        print(f"Plots saved to {savepath}")
        plt.close()
    else:
        plt.show()

def plot_timeseries_model_list(model_list,  start, end, numpoints, savepath):
    if isinstance(model_list[0], str): # If its a list of antimony strings:
        models = []
        for astr in model_list:
            r = te.loada(astr)
            models.append(r)
    else:
        models = model_list
    # Otherwise presumably it's a list of roadrunner models
    n = len(models)
    rows, cols = get_best_dimensions(n)
    idx = 0
    plt.rcParams.update({'font.size': 6})
    fig, axs = plt.subplots(rows, cols)
    for i in range(rows):
        for j in range(cols):
            r = models[idx]
            m = r.simulate(start, end, numpoints)
            axs[i, j].plot(m['time'], m[:, 1:])
            axs[i, j].set_title(idx)
            idx += 1
            if idx > n - 1:
                break
        if idx > n - 1:
            break
    fig.tight_layout()
    if savepath:
        plt.savefig(savepath)
        print(f"Plots saved to {savepath}")
        plt.close()
    else:
        plt.show()

def plot_single_model(model, start, end, numpoints, savepath):
    if isinstance(model, str): # if it's an antimony model
        r = te.loada(model)
    else: #  If it's already a roadrunner model
        r = model
    m = r.simulate(start, end, numpoints)
    r.plot(savefig=savepath)

def plot_timeseries(input, start=0, end=1, numpoints=200, savepath=None):
    if isinstance(input, str) and "->" not in input:  # If we're given a path
        plot_timeseries_path(input, start, end, numpoints, savepath)
    elif isinstance(input, list): # List of models, either as antimony strings, or roadrunner models
        plot_timeseries_model_list(input, start, end, numpoints, savepath)
    else: # antimony string for single model or single roadrunner model
        plot_single_model(input, start, end, numpoints, savepath)

#--------------------------------------------------
# Plotting Fitness Data
#-------------------------------------------------