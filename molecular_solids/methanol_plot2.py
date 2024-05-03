import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def import_arbitrary_module(module_name,path):
    import importlib.util
    import sys

    spec = importlib.util.spec_from_file_location(module_name,path)
    imported_module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = imported_module
    spec.loader.exec_module(imported_module)

    return imported_module

lib         = import_arbitrary_module("loading_lib","/home/users/mendla/work/lib.py")



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


    _start = 53
    _end = 2000
    data = np.load("../scripts/test_2.npy")[_start:_end,:]
    # data = data[~np.ma.fix_invalid(data).mask.any(axis=1)]
    summed = np.flip(np.cumsum(np.flip(data, axis=0),axis=0),axis=0)

    # summed = np.cumsum(data[:,::-1],axis=0)[:,::-1]
    # print(summed)

    # fig, ax = plt.subplots(figsize=(2.79077651389*2,2.6),constrained_layout=True)
    # plt.plot(data[30:,0],summed[30:,2],label="SCF",linewidth=0.5)
    # plt.plot(data[30:,0],summed[30:,2]+summed[30:,3],label="MP2",linewidth=0.5)
    # plt.plot(data[30:,0],summed[30:,1],label="DMA",linewidth=0.5)
    # # plt.yscale("symlog")
    # plt.ylabel("Hartree")
    # plt.xlabel("$a_0$")
    # plt.legend()
    # plt.savefig("CO2_sum.png",dpi=300)

    lib.list_table_data(data,summed,0,_end-_start-1)
    

    fig, ax = plt.subplots(figsize=(2.79077651389*2,5),constrained_layout=True)
    plt.plot(data[:,0],np.abs(summed[:,1]-summed[:,2]),label=r"up to $0^\mathrm{th}$ order",linewidth=0.5)
    plt.plot(data[:,0],np.abs(summed[:,1]-summed[:,3]),label=r"up to $1^\mathrm{st}$ order",linewidth=0.5)
    plt.plot(data[:,0],np.abs(summed[:,1]-summed[:,4]),label=r"up to $2^\mathrm{nd}$ order",linewidth=0.5)
    plt.plot(data[:,0],np.abs(summed[:,1]-summed[:,5]),label=r"up to $3^\mathrm{rd}$ order",linewidth=0.5)
    plt.plot(data[:,0],np.abs(summed[:,1]-summed[:,6]),label=r"up to $4^\mathrm{th}$ order",linewidth=0.5)
    plt.plot(data[:,0],np.abs(summed[:,1]-summed[:,7]),label=r"up to $5^\mathrm{th}$ order",linewidth=0.5)

    # plt.plot(data[30:,0],np.abs(summed[30:,3]),label=r"$|\mathrm{SCF}-\mathrm{MP2}|$",linewidth=0.5)
    plt.yscale("log")
    plt.ylabel("Two body energy / Hartree")
    plt.xlabel("Cutoff distance / $a_0$")
    plt.legend()
    plt.savefig("methanol_orders_final.png",dpi=300)



    # fig, ax = plt.subplots(figsize=(2.79077651389*2,2.6),constrained_layout=True)
    # plt.plot(data[30:,0],np.abs(data[30:,1] - data[30:,2]),label=r"$|\mathrm{SCF}-\mathrm{DMA}|$",linewidth=0.5)
    # plt.plot(data[30:,0],np.abs(data[30:,3]),label=r"$|\mathrm{SCF}-\mathrm{MP2}|$",linewidth=0.5)
    # plt.yscale("symlog")
    # plt.ylabel("Hartree")
    # plt.xlabel("$a_0$")
    # plt.legend()
    # plt.savefig("CO2_differences.png",dpi=300)







    
