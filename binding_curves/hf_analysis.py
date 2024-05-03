#! /home/users/mendla/custom_software/env/bin/python


import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy.optimize import curve_fit

from matplotlib import ticker

fGEOMETRY_1 = np.array([
    [0.0,0.0,0.0],
    [0.917,0.0,0.0]
])
fGEOMETRY_2 = np.array([
    [1.834,0.0,0.0],
    [2.751,0.0,0.0]
])

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

def process_multipoles(data1,data2,i):
    vector,rotation = geometry.LineScanHF01.transform(i)

    for j,entry in enumerate(data2):
        entry.position = (fGEOMETRY_2[j] + vector)/0.529177
    print(vector)

    multipole0 = lib.DMA.energy_between_gdma(data1,data2,0)
    multipole1 = lib.DMA.energy_between_gdma(data1,data2,1) + multipole0
    multipole2 = lib.DMA.energy_between_gdma(data1,data2,2) + multipole1
    multipole3 = lib.DMA.energy_between_gdma(data1,data2,3) + multipole2
    multipole4 = lib.DMA.energy_between_gdma(data1,data2,4) + multipole3
    multipole5 = lib.DMA.energy_between_gdma(data1,data2,5) + multipole4

    return [multipole0,multipole1,multipole2,multipole3,multipole4,multipole5]



def process_file(path,i):
    print(f"Porcessing file \'{path}\'")

    sys.stdout.flush()

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


    data1 = lib.DMA.get_from_gdam_file("../HF_model_data/HF_gdma/0/dma4.punch")
    data2 = lib.DMA.get_from_gdam_file("../HF_model_data/HF_gdma/0/dma4.punch")

    for i,entry in enumerate(data1):
        entry.position = fGEOMETRY_1[i]/0.529177 
    for i,entry in enumerate(data2):
        entry.position = fGEOMETRY_2[i]/0.529177 


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
    for i in geometry.LineScanHF01.generator():
        x.append(geometry.LineScanCO201.plot_position(i))
        CC1,HF1,energy_from_monomers1,_ = process_file(f"../HF2/{format(i+1,'02')}_AVQZ/molpro.out",i)
        CC2,HF2,energy_from_monomers2,_ = process_file(f"../HF2/{format(i+1,'02')}_AV5Z/molpro.out",i)

        multipoles = process_multipoles(data1,data2,i)



        # sapt_electrostatic_energy.append(process_sapt_calculation(f"../HF/{format(i+1,'02')}_AV5Z/molpro.out"))
        
        for array,value in zip(M,multipoles):
            array.append(value)
        
        E.append(lib.get_parameters(5,[CC1[0],CC2[0]])[0])
        E_BSSE.append(lib.get_parameters(5,[CC1[1],CC2[1]])[0])
        E_HF.append(lib.get_parameters(5,[HF1[0],HF2[0]])[0])
        E_HF_BSSE.append(lib.get_parameters(5,[HF1[1],HF2[1]])[0])

    lines_=[]
    _from = 4
    fig,ax = plt.subplots(figsize=(2.79077651389*1.5,2*2*0.9*0.85),constrained_layout=True)


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
        
    plt.savefig(f"HF_GDMA.png",dpi=300)


    with lib.PlottingContextManager(f"HF",[".png",".pdf"]) as (fig,ax):
        for i,values in enumerate(M):
            print(i,values)
            lib.plot_serie(x[_from:],np.array(values[_from:]),label=f"multipoles up to order {i}",color=None)
        # plt.scatter(x[_from:],E_FROM_MON[_from:],label="monomers (old)")
        plt.xlabel("Angstrom")
        lib.plot_serie(x[_from:],E_HF[_from:],label="HF",color=None)
        lib.plot_serie(x[_from:],E[_from:],label="CCSD(T)",color=None)
    
 




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



    





