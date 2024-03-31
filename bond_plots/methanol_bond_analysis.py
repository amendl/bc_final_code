#! /home/users/mendla/custom_software/env/bin/python


import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy.optimize import curve_fit
import seaborn as sns


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



def process_file(path,i,data1,data2):
    multipole = []
    vector,rotation = geometry.LineScanOHBond01.transform(i)

    for j,entry in enumerate(data2):
        entry.position = (geometry.GEOMETRY_3[j] + vector)/0.529177

    
    print(f"Porcessing file \'{path}\'")

    sys.stdout.flush()

    vector,rotation = geometry.LineScanOHBond01.transform(i)
    transformed_geometry_3 = geometry.move_molecule(geometry.rotate_molecule(geometry.GEOMETRY_3,rotation,1),vector)

    lines = lib.read_lines(path)
    
    CC_CP = lib.get_value(lines,"DE_CC_CPAB")
    CC    = lib.get_value(lines,"DE_CC_AB")

    HF = lib.get_values(lines,"!RHF STATE 1.1 Energy","   ")
    
    
    multipole0 = lib.DMA.energy_between_gdma(data1,data2,0)
    multipole1 = lib.DMA.energy_between_gdma(data1,data2,1) + multipole0
    multipole2 = lib.DMA.energy_between_gdma(data1,data2,2) + multipole1
    multipole3 = lib.DMA.energy_between_gdma(data1,data2,3) + multipole2
    multipole4 = lib.DMA.energy_between_gdma(data1,data2,4) + multipole3
    multipole5 = lib.DMA.energy_between_gdma(data1,data2,5) + multipole4


    multipole  = [multipole0,multipole1,multipole2,multipole3,multipole4,multipole5]




    return (CC_CP,CC),(HF[0]-HF[1]-HF[2],HF[0]-HF[3]-HF[4]),multipole



if __name__=="__main__":
    sns.set_theme("paper","ticks")

    params = {#'text.usetex' : True,
            'font.size' : 11.74983
            #   'font.family' : 'lmodern',
            #   'text.latex.unicode': True,
            }
    plt.rcParams.update(params)
    plt.rcParams["font.family"] = "Nimbus Roman"

    # CC2,HF2,c2,d2,e2,f2,energy_from_monomers2,multipoles =  process_file('../oh_bond_2_dipole/01_AVTZ/molpro.out',0)
    # print(multipoles[0])
    # print(multipoles[1])
    # exit()
    x = []
    E = []
    E_BSSE = []
    E_HF = []
    E_HF_BSSE = []

    data1 = lib.DMA.get_from_gdam_file("../methanol_model_data/HF_gdma_V5Z/0/dma4.punch")
    for i,entry in enumerate(data1):
        entry.position = geometry.GEOMETRY_1[i]/0.529177 
    data3 = lib.DMA.get_from_gdam_file("../methanol_model_data/HF_gdma_V5Z/2/dma4.punch")
    for i,entry in enumerate(data3):
        entry.position = geometry.GEOMETRY_3[i]/0.529177 

    M = [[] for _ in range(6)]
    for i in geometry.LineScanOHBond01.generator():
        x.append(geometry.LineScanOHBond01.plot_position(i))
        CC1,HF1,multipoles = process_file(f"../oh_bond_2_dipole/{format(i+1,'02')}_AVTZ/molpro.out",i,data1,data3)
        CC2,HF2,_ = process_file(f"../oh_bond_2_dipole/{format(i+1,'02')}_AVQZ/molpro.out",i,data1,data3)

        for array,value in zip(M,multipoles):
            array.append(value)
        
        E.append(lib.get_parameters(4,[CC1[0],CC2[0]])[0])
        E_BSSE.append(lib.get_parameters(4,[CC1[1],CC2[1]])[0])
        E_HF.append(lib.get_parameters(4,[HF1[0],HF2[0]])[0])
        E_HF_BSSE.append(lib.get_parameters(4,[HF1[1],HF2[1]])[0])
        
    _from = 4
    # with lib.PlotingContextManager("orders_methanol",[".png",".pdf"]) as (fig,ax):
    fig,ax = plt.subplots(figsize=(2.79077651389*1.5,2*2*0.9),constrained_layout=True)

    lines_ = []
    for i,values in enumerate(M):
        print(i,values)
        lines_.append(lib.plot_serie(x[_from:],np.array(values[_from:]),label=f"dma order {i}",color=helper.bc_colors[i]))
    # plt.scatter(x[_from:],E_FROM_MON[_from:],label="monomers (old)")
    # lib.plot_serie(x[_from:],sapt_electrostatic_energy[_from:],label="Electrostatic sapt",color=None)
    plt.xlabel(r"$\AA$")
    plt.ylabel("Hartree")
    lines_.append(lib.plot_serie(x[_from:],E_HF[_from:],label="HF",color=helper.bc_colors["SCF"]))
    lines_.append(lib.plot_serie(x[_from:],E[_from:],label="CCSD(T)",color=helper.bc_colors["CCSD(T)"]))

    plt.savefig("orders_methanol.pdf")
    plt.clf()

    fig,ax = plt.subplots(figsize=(2.79077651389/2,2*2*0.9),constrained_layout=True)
    fig.legend(lines_, [f"DMA order {i}" for i in range(6)] + ["SCF","CCSD(T)"],loc="center")
    plt.axis('off')
    plt.savefig("orders_methanol_legend.pdf")
    
    sys.stdout.flush()
  

