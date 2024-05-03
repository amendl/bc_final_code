#! /home/users/mendla/custom_software/env/bin/python


import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import math

from scipy.optimize import curve_fit
import seaborn as sns
from spherical.wigner import quaternionic


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

def into_real(l,s):
    return np.real([math.sqrt(0.5)*(-s[l+1]+s[l-1]),s[l],math.sqrt(0.5)*(-s[l+1]-s[l-1])/1j])

def get_DMA_from_molpro(path):
    lines = lib.read_lines(path)

    for index in lib.DMA.find(lines,1,path):
        return lib.DMA.get(lines,index,path)
    raise ValueError("Error")
    

if __name__=="__main__":
    sns.set_theme("paper","ticks")
    # #Direct input 
    # plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
    # #Options
    # params = {'text.usetex' : True,
    #       'font.size' : 11,
    #       'font.family' : 'lmodern'
    #       # 'text.latex.unicode': True,
    #       }
    # plt.rcParams.update(params) 


    params = {#'text.usetex' : True,
            'font.size' : 11.74983
            #   'font.family' : 'lmodern',
            #   'text.latex.unicode': True,
            }
    plt.rcParams.update(params)
    plt.rcParams["font.family"] = "Nimbus Roman"

    parser = argparse.ArgumentParser(prog="",description="",epilog="")
    parser.add_argument('-f','--folder',required=True,help='')
    parser.add_argument('-o','--order',required=True,help='')

    args = parser.parse_args()
    folder = args.folder
    order = int(args.order)

    x = np.zeros(90)
    y = np.zeros(90)

    
    with open(f"{folder}/angles") as file:
        for idx in range(30):
            angles = [float(i) for i in file.readline().split(",")]
            DMA_entry = get_DMA_from_molpro(f"{folder}/{idx}/molpro.out")[0]
            DMA_entry_original = get_DMA_from_molpro(f"{folder}/original/molpro.out")[0]
            DMA_entry_original.rotate(quaternionic.array.from_rotation_matrix(lib.create_rotation_matrix(angles[0],angles[1],angles[2])))
            orig = into_real(order,DMA_entry_original.data[order])
            x[idx*3] = orig[0]
            x[idx*3+1] = orig[0+1]
            x[idx*3+2] = orig[0+2]
            new = into_real(order,DMA_entry.data[order])
            y[idx*3] = new[0]
            y[idx*3+1] = new[0+1]
            y[idx*3+2] = new[0+2]

        fig,ax = plt.subplots(figsize=(2.79077651389,2*2*0.9),constrained_layout=True)

        plt.scatter(x,y,s=5,label="Coefficients",zorder=1)
        plt.plot(np.linspace(np.min(x),np.max(x),100),np.linspace(np.min(y),np.max(y),100),linestyle="dashed",color="red",label="Ideal",zorder=0)

        if order==1:
            plt.xlabel(r"Rotated molecule / $ea_0$")
            plt.ylabel(r"Rotated spherical expansion / $ea_0$")
        else:
            plt.xlabel(f'Rotated molecule / $ea_0^{order}$')
            plt.ylabel(f'Rotated spherical expansion / $ea_0^{order}$')

        plt.legend()



            

        plt.savefig(f"molpro_{folder}_{order}.png",dpi=300)


            
    


