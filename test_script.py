import tellurium as te

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os
import pandas as pd

import math

def plot_several_models_and_savefig(filepath, outpath, rows, cols):
    total = len(os.listdir(filepath))
    num_rounds = math.ceil(total/(rows*cols))
    idx = 0
    models = []
    model_names = []
    for dir in os.listdir(filepath):
        parentdir = os.path.join(filepath, dir)
        fullpath = os.path.join(parentdir, "final_models")
        for model_file in os.listdir(fullpath):
            if model_file.startswith("bestmodel"):
                f = open(os.path.join(fullpath, model_file), "r")
                astr = f.read()
                f.close()
                models.append(astr)
                model_names.append(parentdir[-15:])
    for rnd in range(num_rounds):

        fig, axs = plt.subplots(rows, cols)
        for i in range(rows):
            if idx >= total:
                plt.savefig(os.path.join(outpath, f"{rnd}"), dpi=5000)
                plt.clf()
                break
            for j in range(cols):
                r = te.loada(models[idx])
                try:
                    m = r.simulate(0)
                    r.steadyState()
                    axs[i, j].plot(m['time'], m[:, 1:])
                    axs[i, j].set_title(f"{model_names[idx]}\n{has_oscillator_eigens(r.getFullEigenValues())}", size=8)
                except:
                    pass
                idx += 1
                if idx >= total:
                    plt.savefig(os.path.join(outpath, f"{rnd}"), dpi=5000)
                    plt.clf()
                    break
        plt.show()
        plt.savefig(os.path.join(outpath, f"{rnd}"))
        plt.clf()
    print("done")



def plot_several_models(filepath, rows, cols):
    fig, axs = plt.subplots(rows, cols)
    models = []
    model_names = []
    cwd = os.getcwd()
    cwd = os.path.join(cwd, "src")
    for file in os.listdir(filepath):
        model_names.append(file[6:])
        fullpath = os.path.join(filepath, file)
        fullpath = os.path.join(cwd, fullpath)
        f = open(fullpath, "r")
        models.append(f.read())
        f.close()
    idx = 0
    for i in range(rows):
        if idx >= len(model_names):
            break
        for j in range(cols):
            r = te.loada(models[idx])
            try:
                m = r.simulate(0)
                r.steadyState()
                axs[i,j].plot(m['time'], m[:,1:])
                axs[i,j].set_title(f"{model_names[idx]} {has_oscillator_eigens(r.getFullEigenValues())}", size=8)
            except:
                pass
            idx += 1
            if idx >= len(models):
                break
    plt.show()

def has_oscillator_eigens(eigenvals):
    for num in eigenvals:
        if num.real > 0 and num.imag != 0:
            return True
    return False


# path ="/home/hellsbells/Desktop/networkEv/src/final_models"
# plot_several_models(path, 5, 2)

def make_truth_csv(model, filename, simulation=[0, 1, 10]):
    if not filename.endswith(".csv"):
        filename += ".csv"
    cwd = os.getcwd()
    path = os.path.join(cwd, f"test_files/{filename}")
    r = te.loada(model)
    result = r.simulate(simulation[0], simulation[1], simulation[2])
    df = pd.DataFrame(result, columns=['time']+r.getFloatingSpeciesIds())
    df.set_index('time', inplace=True)
    df.to_csv(path)

def plot_single_model(model, simulation=[0, 1, 10]):
    r = te.loada(model)
    r.simulate(simulation[0], simulation[1], simulation[2])
    r.plot()

def compare_model_plots(truthmodel, evolvedmodel, simulation=[0, 1, 10]):
    r1 = te.loada(truthmodel)
    result1 = r1.simulate(simulation[0], simulation[1], simulation[2])

    r2 = te.loada(evolvedmodel)
    result2 = r2.simulate(simulation[0], simulation[1], simulation[2])

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1, title="Truth Model")
    ax2 = fig.add_subplot(1, 2, 2, title="Evolved Model")#, sharey=ax1)
    ax2.set_ylim(0,20)
    ax1.plot(result1[:, 0], result1[:, 1:])
    ax2.plot(result2[:, 0], result2[:, 1:])
    ax1.legend(r1.getFloatingSpeciesIds())
    ax2.legend(r2.getFloatingSpeciesIds())
    plt.show()


model = """S2 -> S1; k0*S2
S2 + S0 -> S0; k1*S2*S0
S0 -> S2; k2*S0
S2 -> S2+S2; k3*S2
S1 -> S0+S2; k4*S1
k0 = 3.708920735583032
k1 = 4.027285355118902
k2 = 7.450108022248875
k3 = 145.89432306715088
k4 = 16.371587764670334
S0 = 1.0
S1 = 5.0
S2 = 9.0

"""

evolvedmodel = """
// Created by libAntimony v2.12.0
// Compartments and Species:
species S0, S1, S2;

// Reactions:
_J0: S0 + S1 -> S1; k1*S0*S1;
_J1: S1 + S2 -> S1 + S2; k2*S1*S2;
_J2: S0 + S2 -> S1 + S1; k3*S0*S2;
_J3: S0 -> S1 + S1; k4*S0;
_J4: S0 -> S1; k5*S0;
_J5: S1 + S2 -> S1 + S1; k6*S1*S2;
_J6: S0 + S2 -> S0 + S1; k7*S0*S2;
_J7: S2 -> S1 + S2; k8*S2;
_J8: S2 -> S1 + S1; k9*S2;
_J9: S2 -> S2 + S2; k10*S2;
_J10: S2 + S2 -> S1; k11*S2*S2;
_J11: S2 + S2 -> S0; k12*S2*S2;
_J12: S2 -> S0 + S1; k13*S2;
_J13: S1 + S1 -> S1 + S1; k14*S1*S1;
_J14: S0 + S1 -> S0; k15*S0*S1;
_J15: S2 + S2 -> S0 + S2; k16*S2*S2;
_J16: S0 -> S0 + S1; k17*S0;
_J17: S2 -> S0 + S2; k18*S2;
_J18: S1 -> S0 + S1; k19*S1;
_J19: S2 -> S2; k20*S2;
_J20: S0 + S2 -> S1; k21*S0*S2;
_J21: S2 -> S1; k22*S2;
_J22: S2 -> S0; k23*S2;
_J23: S1 + S2 -> S1; k24*S1*S2;

// Species initializations:
S0 = 1;
S1 = 5;
S2 = 9;

// Variable initializations:
k1 = 6.41994618410703;
k2 = 0.0783887558898663;
k3 = 11.6228695729583;
k4 = 4.71537624468825;
k5 = 12.2514262356933;
k6 = 13.6248212876357;
k7 = 30.1783681699379;
k8 = 27.1515698331225;
k9 = 7.98561518156907;
k10 = 3.06230146453495;
k11 = 4.79268604511165;
k12 = 7.72034888874054;
k13 = 11.229991807324;
k14 = 7.36319961101907;
k15 = 7.30870325529001;
k16 = 7.99995231066672;
k17 = 14.2199701328281;
k18 = 4.61755712855557;
k19 = 3.33735072232818;
k20 = 4.29231822957825;
k21 = 1.17627121315487;
k22 = 3.52364195550499;
k23 = 6.01086830929271;
k24 = 1.78126089564575;
"""

# compare_model_plots(model, evolvedmodel, simulation=[0, 10, 1000])
#
# r = te.loada(evolvedmodel)
# r.simulate(0,1000,5000)
# r.plot()
# r.steadyState()
# print(r.getFullEigenValues())
# # print(r.getFullEigenValues())

input= "/home/hellsbells/Desktop/networkEv/Data/BenchmarkTests/1000_trials"
output = "/home/hellsbells/Desktop/networkEv/Data/BenchmarkTests/1000_trials_oscillators"
plot_several_models_and_savefig(input, output, 5, 5)

# plot_several_models("/home/hellsbells/Desktop/networkEv/Data/BenchmarkTests/1000_trials_oscillators",2, 11)


# path = "/home/hellsbells/Desktop/networkEv/Data/BenchmarkTests/1000_trials_oscillators"
# modelname = "bestmodel_H4Pyg4D9SZ5i"
#
# model_path = os.path.join(path, modelname)
#
# f = open(model_path, "r")
# astr = f.read()
# f.close()
# r = te.loada(astr)
# r.simulate(0, 1, 100)
# print(r.getFullEigenValues())
# r.plot()

r = te.loada("""
// Created by libAntimony v2.12.0
// Compartments and Species:
species S0, S2, S1;

// Reactions:
_J0: S0 -> S0 + S0; k1*S0;
_J1: S0 + S2 -> S1 + S1; k2*S0*S2;
_J2: S2 -> S1 + S2; k3*S2;
_J3: S2 -> S2 + S2; k4*S2;
_J4: S2 -> S1 + S1; k5*S2;
_J5: S0 -> S0 + S2; k6*S0;
_J6: S0 + S1 -> S0; k7*S0*S1;
_J7: S1 -> S0 + S1; k8*S1;
_J8: S2 -> S0 + S2; k9*S2;
_J9: S1 -> S1 + S1; k10*S1;
_J10: S0 + S2 -> S1; k11*S0*S2;

// Species initializations:
S0 = 1;
S2 = 9;
S1 = 2;

// Variable initializations:
k1 = 4.21910131567945;
k2 = 11.712073055039;
k3 = 11.0259246184376;
k4 = 8.88687795278714;
k5 = 2.40725961099904;
k6 = 7.46683140994304;
k7 = 2.91322140978341;
k8 = 5.68283591047148;
k9 = 7.63267449265815;
k10 = 13.0555066455122;
k11 = 8.77611525023671;

// Other declarations:
const k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11;""")
# r.simulate()
# r.plot()
#
# r.steadyState()
# print(r.getFullEigenValues())
# r.reset()
# fig, axs = plt.subplots(1, 3)
#
#
# print(r.getFullEigenValues())
# m = r.simulate(0,10,500)
# axs[0].plot(m['time'], m[:,1:])
# axs[0].set_title(f"{r.getFullEigenValues()}", size=8)
#
# r.reset()
# m = r.simulate(0,1000,5000)
# axs[1].plot(m['time'], m[:,1:])
# axs[1].set_title(f"{r.getFullEigenValues()}", size=8)
#
# r.reset()
# m = r.simulate(0,10000,50000)
# axs[2].plot(m['time'], m[:,1:])
# axs[2].set_title(f"{r.getFullEigenValues()}", size=8)
# plt.show()