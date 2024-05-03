#! /home/users/mendla/custom_software/env/bin/python


import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns


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

	
    # POP_values = []
    # for i,index in enumerate(lib.POP.find(lines,5,path)):
    #     POP_values.append(lib.POP.get(lines,index,6 if i in [0,1,2] else 3,path))
        
    
    DMA_entries = []
    for index in lib.DMA.find(lines,5,path):
        DMA_entries.append(lib.DMA.get(lines,index,path))
    
    # print(DMA_entries[0][0])
    
    # E0 = lib.POP.energy(np.concatenate((geometry.GEOMETRY_1,transformed_geometry_3)),POP_values[0])
    # E1 = lib.POP.energy(np.concatenate((geometry.GEOMETRY_1,transformed_geometry_3)),POP_values[1])
    # E2 = lib.POP.energy(np.concatenate((geometry.GEOMETRY_1,transformed_geometry_3)),POP_values[2])
    # E3 = lib.POP.energy(geometry.GEOMETRY_1,POP_values[3])
    # E4 = lib.POP.energy(transformed_geometry_3,POP_values[4])

    energy_from_monomers = lib.DMA.energy_old(DMA_entries[3]+DMA_entries[4],[(len(DMA_entries[3])-1,len(DMA_entries[3]) + len(DMA_entries[4])-1)]) - lib.DMA.energy_old(DMA_entries[3]) - lib.DMA.energy_old(DMA_entries[4]) 
    # multipole = [lib.DMA.energy(DMA_entries[3]+DMA_entries[4],[(len(DMA_entries[3])-1,len(DMA_entries[3]) + len(DMA_entries[4])-1)],[],i) - lib.DMA.energy(DMA_entries[3],[],[],i) - lib.DMA.energy(DMA_entries[4],[],[],i) for i in range(1,3)]
    
    # multipole0 = lib.DMA.energy_between(DMA_entries[3],DMA_entries[4],0,[len(DMA_entries[3])-1],[len(DMA_entries[4])-1])
    # multipole1 = lib.DMA.energy_between(DMA_entries[3],DMA_entries[4],1,[len(DMA_entries[3])-1],[len(DMA_entries[4])-1]) + multipole0
    # multipole2 = lib.DMA.energy_between(DMA_entries[3],DMA_entries[4],2,[len(DMA_entries[3])-1],[len(DMA_entries[4])-1]) + multipole1
    # multipole3 = lib.DMA.energy_between(DMA_entries[3],DMA_entries[4],3,[len(DMA_entries[3])-1],[len(DMA_entries[4])-1]) + multipole2

    multipole = []
    last = 0.
    for i in range(6):
        multipole.append(lib.DMA.energy_between([DMA_entries[0][0],DMA_entries[0][1]],[DMA_entries[0][2],DMA_entries[0][3]],i)+last)
        last = multipole[-1]
        
    # multipole0 = lib.DMA.energy_between([DMA_entries[0][0],DMA_entries[0][1]],[DMA_entries[0][2],DMA_entries[0][3]],0)
    # multipole1 = lib.DMA.energy_between([DMA_entries[0][0],DMA_entries[0][1]],[DMA_entries[0][2],DMA_entries[0][3]],1) + multipole0
    # multipole2 = lib.DMA.energy_between([DMA_entries[0][0],DMA_entries[0][1]],[DMA_entries[0][2],DMA_entries[0][3]],2) + multipole1
    # multipole3 = lib.DMA.energy_between([DMA_entries[0][0],DMA_entries[0][1]],[DMA_entries[0][2],DMA_entries[0][3]],3) + multipole2


    dipole_lines = lib.find_in_lines(lines, "!RHF STATE 1.1 Dipole moment")
    # print(lib.DMA.dipole(DMA_entries[3]))
    # print(lines[dipole_lines[3]])
    # print(lines[dipole_lines[3]+1])

    # print(lib.DMA.dipole(DMA_entries[4]))
    # print(lines[dipole_lines[4]])
    # print(lines[dipole_lines[4]+1])


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

    # CC2,HF2,c2,d2,e2,f2,energy_from_monomers2,multipoles =  process_file('../oh_bond_2_dipole/01_AVTZ/molpro.out',0)
    # print(multipoles[0])
    # print(multipoles[1])
    # exit()
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
        # lib.plot_serie(x[_from:],E_HF[_from:],label="HF",color="black")
        for i,values in enumerate(M):
            print(i,values)
            lines_.append(lib.plot_serie(np.array(x[_from:])+1.74,np.array(values[_from:]),label=f"DMA order {i}",color=helper.bc_colors[i]))
        # plt.scatter(np.array(x[_from:]),E_FROM_MON[_from:],label="monomers (old)")
        plt.xlabel(r"Distance / $\mathrm{\AA}$")
        plt.ylabel("Intermolecular energy / Hartree")
        lines_.append(lib.plot_serie(np.array(x[_from:])+1.74,E_HF[_from:],label="HF",color=helper.bc_colors["SCF"]))
        lines_.append(lib.plot_serie(np.array(x[_from:])+1.74,E[_from:],label="CCSD(T)",color=helper.bc_colors["CCSD(T)"]))

        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1,1))
        ax.yaxis.set_major_formatter(formatter)
        
        plt.savefig(f"orders_HF_{c2}V{c}Z.png",dpi=300)
        plt.clf()
        plt.cla()


    print(lines_)
    fig,ax = plt.subplots(figsize=(2.79077651389*2,0.9),constrained_layout=True)
    fig.legend(lines_, [f"DMA order {i}" for i in range(6)] + ["SCF","CCSD(T)"],loc="center",ncol=3)
    plt.axis('off')
    plt.savefig("orders_HF_legend.png",dpi=300)

        

    # with lib.PlottingContextManager("HF_comparison.pdf") as  (fig, ax):
    #     lib.plot_serie(x,E_HF,label="E HF",color='orange')
    #     lib.plot_serie(x,E_HF_BSSE,label="E HF BSSE",color='yellow')
    #     lib.plot_serie(x,P,label="POP",color="blue")
    #     # lib.plot_serie(x,P_BSSE,label="POP BSSE",color="red")
    #     lib.plot_serie(x,M,label="DMA",color="black")
    #     lib.plot_serie(x,M_BSSE,label="DMA BSSE",color="grey")
    #     lib.plot_serie(x,E_FROM_MON,label="From monomers",color="red")
    #     # lib.plot_serie(x,A,label="charge",color='blue')

    # with lib.PlottingContextManager("CC_comparison.pdf") as  (fig, ax):
    #     lib.plot_serie(x,E,label="E CC",color='orange')
    #     lib.plot_serie(x,E_BSSE,label="E CC BSSE",color='yellow')
    #     lib.plot_serie(x,P,label="POP",color="blue")
    #     # lib.plot_serie(x,P_BSSE,label="POP BSSE",color="red")
    #     lib.plot_serie(x,M,label="DMA",color="black")
    #     lib.plot_serie(x,M_BSSE,label="DMA BSSE",color="grey")
    #     lib.plot_serie(x,E_FROM_MON,label="From monomers",color="red")

    # with lib.PlottingContextManager("check_DMA.pdf") as (fig,ax):
    #     lib.plot_serie(x,E,label="E CC",color='orange')
    #     lib.plot_serie(x,M,label="DMA",color="black")
    #     lib.plot_serie(x,M_BSSE,label="DMA BSSE",color="grey")

    # with lib.PlottingContextManager("monomers_comparison",[".png",".pdf"]) as (fig,ax):
    #     lib.plot_serie(x,E_HF,label="E HF",color="grey")
    #     lib.plot_serie(x,E_HF_BSSE,label="E HF BSSE",color="black")
    #     lib.plot_serie(x,E,label="E CCSD(T)",color="cyan")
    #     lib.plot_serie(x,E,label="E CCSD(T)",color="blue")
    #     lib.plot_serie(x,E_FROM_MON,label="From Monomers",color="red")
    
 




    # def func(x,a,b,c):
    #     return a*np.power(x-b,c)
    
    # popt, pcov = curve_fit(func,x,M[0],p0=[1.,0.,-2.])
    # print(popt[2])
    # popt, pcov = curve_fit(func,x,M[1],p0=[1.,0.,-2.])
    # print(popt[2])
    # popt, pcov = curve_fit(func,x,M[2],p0=[1.,0.,-2.])
    # print(popt[2])
    # popt, pcov = curve_fit(func,x,M[3],p0=[1.,0.,-2.])
    # print(popt[2])

    sys.stdout.flush()


        

    # print(E[-1])
    # print(P[-1])
    # print(M[-1])
    # print(E_FROM_MON)



    





