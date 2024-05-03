import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

if __name__=="__main__":

    sns.set_theme("paper","ticks")

    params = {#'text.usetex' : True,
            'font.size' : 11.74983
            #   'font.family' : 'lmodern',
            #   'text.latex.unicode': True,
            }
    plt.rcParams.update(params)
    plt.rcParams["font.family"] = "Nimbus Roman"

    np.set_printoptions(threshold=np.inf)



    data = np.load("../scripts/test_Me.npy")[14:34,:]
    # print(data)
    # data = data[~np.ma.fix_invalid(data).mask.any(axis=1)]
    summed = np.flip(np.cumsum(np.flip(data, axis=0),axis=0),axis=0)

    # summed = np.cumsum(data[:,::-1],axis=0)[:,::-1]
    # print(summed)

    # _from = 14
    # _to = 34
    _from = 0
    _to = -1


    fig, ax = plt.subplots(figsize=(2.79077651389,3.8),constrained_layout=True)
    plt.plot(data[_from:_to,0],summed[_from:_to,2],label="SCF",linewidth=0.5)
    # plt.plot(data[_from:_to,0],summed[_from:_to,2]+summed[_from:_to,3],label="MP2",linewidth=0.5)
    plt.plot(data[_from:_to,0],summed[_from:_to,1],label="DMA",linewidth=0.5)
    # plt.yscale("symlog")
    plt.ylabel("Two body energy / Hartree")
    plt.xlabel("Cutoff distance / $a_0$")

    plt.legend()
    plt.savefig("Me_sum_final2.png",dpi=300)

    fig, ax = plt.subplots(figsize=(2.79077651389*2,2.6),constrained_layout=True)
    plt.plot(data[_from:_to,0],np.abs(summed[_from:_to,1]-summed[_from:_to,2]),label=r"$|\mathrm{SCF}-\mathrm{DMA}|$",linewidth=0.5)
    plt.plot(data[_from:_to,0],np.abs(summed[_from:_to,3]),label=r"$|\mathrm{MP2\ correction}|$",linewidth=0.5)
    plt.yscale("log")
    plt.ylabel("Two body energy / Hartree")
    plt.xlabel("Cutoff distance / $a_0$")
    plt.legend()
    plt.savefig("NH3_sum_differences_final.png",dpi=300)

    fig, ax = plt.subplots(figsize=(2.79077651389*2,2.6),constrained_layout=True)
    plt.plot(data[_from:_to,0],np.abs(data[_from:_to,1] - data[_from:_to,2]),label=r"$|\mathrm{SCF}-\mathrm{DMA}|$",linewidth=0.5)
    plt.plot(data[_from:_to,0],np.abs(data[_from:_to,3]),label=r"$|\mathrm{MP2\ correction}|$",linewidth=0.5)
    plt.yscale("log")
    plt.ylabel("Two body energy / Hartree")
    plt.xlabel("Cutoff distance / $a_0$")
    plt.legend()
    plt.savefig("me_differences_final.png",dpi=300)
