import seaborn as sns
import matplotlib.pyplot as plt

import pickle
import numpy as np

def zscores(std, mean, score):


    return (std * score) + mean


def plot_dist(data, data_point, name,label="RMSD (Å)"):

    # fig = plt.figure(figsize=(3,6))

    data = np.array(data)
    std = np.std(data)
    mean = np.average(data)

    sns.displot(data, binwidth=.1,)#kind="kde", cut=0,bw_adjust=.25)
    plt.xlabel(f"{label} \n Gray lines represent -1, 0 (mean) and 1 Z scores\n Red line represents Z score for hit")
    
    plt.axvline(x=zscores(std, mean, -1), color='gray')
    plt.axvline(x=mean, color='gray')
    plt.axvline(x=zscores(std, mean, 1), color='gray')
    if data_point:
        plt.axvline(zscores(std, mean, data_point), color='r')


    fig = plt.gcf()
    fig.set_size_inches( 6, 3)
    plt.tight_layout()
    plt.savefig(name)

    plt.close("all")

def plot_bivariate(lens,rmsds):
    sns.displot(x=lens, y=rmsds, binwidth=(1, .1))
    plt.savefig("plot_dist.png")
    plt.close()
    



if __name__ == "__main__":
    with open("/app/output/EPI_PDB_fragment_pairs_6XR8_A_exp.pickle", 'rb') as inhandle: 
        data  = pickle.load(inhandle)
    plot_bivariate(data["lens"], data["rmsds"])
