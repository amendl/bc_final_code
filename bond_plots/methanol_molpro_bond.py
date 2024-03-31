#! /home/users/mendla/custom_software/env/bin/python


import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy.optimize import curve_fit
import seaborn as sns
from matplotlib import ticker



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

	
    POP_values = []
    for i,index in enumerate(lib.POP.find(lines,5,path)):
        POP_values.append(lib.POP.get(lines,index,12 if i in [0,1,2] else 6,path))
        
    
    DMA_entries = []
    for index in lib.DMA.find(lines,5,path):
        DMA_entries.append(lib.DMA.get(lines,index,path))
    
    # print(DMA_entries[0][0])
    
    E0 = lib.POP.energy(np.concatenate((geometry.GEOMETRY_1,transformed_geometry_3)),POP_values[0])
    E1 = lib.POP.energy(np.concatenate((geometry.GEOMETRY_1,transformed_geometry_3)),POP_values[1])
    E2 = lib.POP.energy(np.concatenate((geometry.GEOMETRY_1,transformed_geometry_3)),POP_values[2])
    E3 = lib.POP.energy(geometry.GEOMETRY_1,POP_values[3])
    E4 = lib.POP.energy(transformed_geometry_3,POP_values[4])

    energy_from_monomers = lib.DMA.energy_old(DMA_entries[3]+DMA_entries[4],[(len(DMA_entries[3])-1,len(DMA_entries[3]) + len(DMA_entries[4])-1)]) - lib.DMA.energy_old(DMA_entries[3]) - lib.DMA.energy_old(DMA_entries[4]) 
    # multipole = [lib.DMA.energy(DMA_entries[3]+DMA_entries[4],[(len(DMA_entries[3])-1,len(DMA_entries[3]) + len(DMA_entries[4])-1)],[],i) - lib.DMA.energy(DMA_entries[3],[],[],i) - lib.DMA.energy(DMA_entries[4],[],[],i) for i in range(1,3)]
    
    multipole0 = lib.DMA.energy_between(DMA_entries[3],DMA_entries[4],0,[len(DMA_entries[3])-1],[len(DMA_entries[4])-1])
    multipole1 = lib.DMA.energy_between(DMA_entries[3],DMA_entries[4],1,[len(DMA_entries[3])-1],[len(DMA_entries[4])-1]) + multipole0
    multipole2 = lib.DMA.energy_between(DMA_entries[3],DMA_entries[4],2,[len(DMA_entries[3])-1],[len(DMA_entries[4])-1]) + multipole1
    multipole3 = lib.DMA.energy_between(DMA_entries[3],DMA_entries[4],3,[len(DMA_entries[3])-1],[len(DMA_entries[4])-1]) + multipole2

    multipole  = [multipole0,multipole1,multipole2,multipole3]

    dipole_lines = lib.find_in_lines(lines, "!RHF STATE 1.1 Dipole moment")
   

    return (CC_CP,CC),(HF[0]-HF[1]-HF[2],HF[0]-HF[3]-HF[4]),E0-E1-E2, E0 - E3-E4,None,None, energy_from_monomers,multipole



if __name__=="__main__":
    sns.set_theme("paper","ticks")

    params = {#'text.usetex' : True,
            'font.size' : 11.74983
            #   'font.family' : 'lmodern',
            #   'text.latex.unicode': True,
            }
    plt.rcParams.update(params)
    plt.rcParams["font.family"] = "Nimbus Roman"

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
    M = [[] for _ in range(4)]
    for i in geometry.LineScanOHBond01.generator():
        x.append(geometry.LineScanOHBond01.plot_position(i))
        CC1,HF1,c1,d1,e1,f1,energy_from_monomers1,_ = process_file(f"../oh_bond_2_dipole/{format(i+1,'02')}_AVTZ/molpro.out",i)
        CC2,HF2,c2,d2,e2,f2,energy_from_monomers2,multipoles = process_file(f"../oh_bond_2_dipole/{format(i+1,'02')}_AVQZ/molpro.out",i)

        sapt_electrostatic_energy.append(process_sapt_calculation(f"../sapt/{format(i+1,'02')}_AVQZ/molpro.out"))
        
        for array,value in zip(M,multipoles):
            array.append(value)
        
        E.append(lib.get_parameters(4,[CC1[0],CC2[0]])[0])
        E_BSSE.append(lib.get_parameters(4,[CC1[1],CC2[1]])[0])
        E_HF.append(lib.get_parameters(4,[HF1[0],HF2[0]])[0])
        E_HF_BSSE.append(lib.get_parameters(4,[HF1[1],HF2[1]])[0])
        P.append(c2)
        P_BSSE.append(d2)
        M_BSSE.append(f2)
        E_FROM_MON.append(energy_from_monomers2)


    _from = 4
    # with lib.PlottingContextManager("new_orders_comparison_from_d_matrix_AVTZ",[".png",".pdf"]) as (fig,ax):
    fig,ax = plt.subplots(figsize=(2.79077651389,2*2*0.9),constrained_layout=True)


    for i,values in enumerate(M):
        print(i,values)
        lib.plot_serie(x[_from:],np.array(values[_from:]),label=f"multipoles up to order {i}",color=helper.bc_colors[i])
    plt.xlabel(r"$\AA$")
    plt.ylabel("Hartree")
    lib.plot_serie(x[_from:],E[_from:],label="CCSD(T)",color=helper.bc_colors["CCSD(T)"])
    lib.plot_serie(x[_from:],E_HF[_from:],label="HF",color=helper.bc_colors["SCF"])


    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1,1))
    ax.yaxis.set_major_formatter(formatter)

    plt.savefig("methanol_bond_molpro_AVQZ.pdf")

    sys.stdout.flush()
    
