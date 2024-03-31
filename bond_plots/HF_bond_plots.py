#! /home/users/mendla/custom_software/env/bin/python


import numpy as np
import matplotlib.pyplot as plt
import sys

from matplotlib import ticker

from scipy.optimize import curve_fit

def import_arbitrary_module(module_name,path):
    import importlib.util
    import sys

    spec = importlib.util.spec_from_file_location(module_name,path)
    imported_module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = imported_module
    spec.loader.exec_module(imported_module)

    return imported_module

lib         = import_arbitrary_module("loading_lib","/home/users/mendla/work/lib.py")
geometry    = import_arbitrary_module("geometries","/home/users/mendla/work/geometries.py")
helper    = import_arbitrary_module("helper","/home/users/mendla/work/helper.py")



def process_sapt_calculation(path):
    print(f"Processing file \'{path}\'")
    sys.stdout.flush()

    lines = lib.read_lines(path)
    energy = lib.get_value(lines,"E_HF")/1000.
    print(energy)
    return energy



def process_file(path,i):
    print(f"Porcessing file \'{path}\'")

    sys.stdout.flush()

    vector,rotation = geometry.LineScanOHBond01.transform(i)
    transformed_geometry_3 = geometry.move_molecule(geometry.rotate_molecule(geometry.GEOMETRY_3,rotation,1),vector)

    lines = lib.read_lines(path)
    
    CC_CP = lib.get_value(lines,"DE_CC_CPAB")
    CC    = lib.get_value(lines,"DE_CC_AB")

    HF = lib.get_values(lines,"!RHF STATE 1.1 Energy","   ")

	
         
    
    DMA_entries = []
    for index in lib.DMA.find(lines,5,path):
        DMA_entries.append(lib.DMA.get(lines,index,path))
    
  
    energy_from_monomers = lib.DMA.energy_old(DMA_entries[3]+DMA_entries[4],[(len(DMA_entries[3])-1,len(DMA_entries[3]) + len(DMA_entries[4])-1)]) - lib.DMA.energy_old(DMA_entries[3]) - lib.DMA.energy_old(DMA_entries[4]) 
 
    multipole = []
    last = 0.
    for i in range(6):
        multipole.append(lib.DMA.energy_between([DMA_entries[0][0],DMA_entries[0][1]],[DMA_entries[0][2],DMA_entries[0][3]],i)+last)
        last = multipole[-1]
     dipole_lines = lib.find_in_lines(lines, "!RHF STATE 1.1 Dipole moment")


    return (CC_CP,CC),(HF[0]-HF[1]-HF[2],HF[0]-HF[3]-HF[4]), energy_from_monomers,multipole



if __name__=="__main__":
    sns.set_theme("paper","ticks")

    params = {#'text.usetex' : True,
            'font.size' : 11.74983
            #   'font.family' : 'lmodern',
            #   'text.latex.unicode': True,
            }
    plt.rcParams.update(params)
    plt.rcParams["font.family"] = "Nimbus Roman"

    c2='A'

    lines_ = []


    for c in ["D","T","Q","5"]:
        x = []
        E = []
        E_BSSE = []
        E_HF = []
        E_HF_BSSE = []
        P = []
        P_BSSE = []
        M = []
        M_BSSE = []
        E_FROM_MON = []
        sapt_electrostatic_energy = []
        M = [[] for _ in range(6)]
        for i in geometry.LineScanCO201.generator():
            x.append(geometry.LineScanCO201.plot_position(i))
            CC1,HF1,energy_from_monomers1,_ = process_file(f"../HF2/{format(i+1,'02')}_AVQZ/molpro.out",i)
            CC2,HF2,energy_from_monomers2,_ = process_file(f"../HF2/{format(i+1,'02')}_AV5Z/molpro.out",i)

            _,_,_,multipoles = process_file(f"../HF2/{format(i+1,'02')}_{c2}V{c}Z/molpro.out",i)


            sapt_electrostatic_energy.append(process_sapt_calculation(f"../HF/{format(i+1,'02')}_AV5Z/molpro.out"))
            
            for array,value in zip(M,multipoles):
                array.append(value)
            
            E.append(lib.get_parameters(5,[CC1[0],CC2[0]])[0])
            E_BSSE.append(lib.get_parameters(5,[CC1[1],CC2[1]])[0])
            E_HF.append(lib.get_parameters(5,[HF1[0],HF2[0]])[0])
            E_HF_BSSE.append(lib.get_parameters(5,[HF1[1],HF2[1]])[0])

        _from = 4
        fig,ax = plt.subplots(figsize=(2.79077651389,2*2*0.9),constrained_layout=True)

        lines_.clear()

        # with lib.PlottingContextManager(f"orders_HF_{c2}V{c}Z",[".png",".pdf"]) as (fig,ax):
        lib.plot_serie(x[_from:],E_HF[_from:],label="HF",color="black")
        for i,values in enumerate(M):
            print(i,values)
            lines_.append(lib.plot_serie(x[_from:],np.array(values[_from:]),label=f"DMA order {i}",color=helper.bc_colors[i]))
        # plt.scatter(x[_from:],E_FROM_MON[_from:],label="monomers (old)")
        plt.xlabel(r"$\AA$")
        plt.ylabel("Hartree")
        lines_.append(lib.plot_serie(x[_from:],E_HF[_from:],label="HF",color=helper.bc_colors["SCF"]))
        lines_.append(lib.plot_serie(x[_from:],E[_from:],label="CCSD(T)",color=helper.bc_colors["CCSD(T)"]))

        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1,1))
        ax.yaxis.set_major_formatter(formatter)
        
        plt.savefig(f"orders_HF_{c2}V{c}Z.pdf")
        plt.clf()
        plt.cla()


    print(lines_)
    fig,ax = plt.subplots(figsize=(2.79077651389*2,0.9),constrained_layout=True)
    fig.legend(lines_, [f"DMA order {i}" for i in range(6)] + ["SCF","CCSD(T)"],loc="center",ncol=3)
    plt.axis('off')
    plt.savefig("orders_HF_legend.pdf")

       

    sys.stdout.flush()
