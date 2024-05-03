#! /home/users/mendla/custom_software/env/bin/python

'''
# Usage:
cd work/2/HF2
/home/users/mendla/work/stability.py --files="01_AVDZ/molpro.out;01_VDZ/molpro.out;01_AVTZ/molpro.out;01_VTZ/molpro.out;01_AVQZ/molpro.out;01_VQZ/molpro.out;01_AV5Z/molpro.out;01_V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="3" --site_index="0" --no_of_records="5" --output="HF_H"
cd work/1/13
/home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z_simple/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="0" --no_of_records="1" --output="methanol_H"
/home/users/mendla/work/stability.py --files="AVDZ/molpro_cc.out;VDZ/molpro_cc.out;AVTZ/molpro_cc.out;VTZ/molpro_cc.out;AVQZ/molpro_cc.out;VQZ/molpro_cc.out;AV5Z_simple/molpro_cc.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="0" --no_of_records="1" --output="methanol_CC"
cd work/2/CO2_molpro_HF
/home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="1" --no_of_records="1" --output="CO2_H"
cd work/2/CO2_molpro
/home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="1" --no_of_records="1" --output="CO2_CC"
cd work/2/HF_molpro_HF
/home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="1" --no_of_records="1" --output="HF_H"
cd work/2/HF_molpro
/home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="1" --no_of_records="1" --output="HF_CC"
cd work/2/NH3_molpro_HF
/home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="0" --no_of_records="1" --output="NH3_H"
cd work/2/NH3_molpro
/home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="0" --no_of_records="1" --output="NH3_CC"
'''

import os
import shutil
import math 
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

# plt.rcParams['text.usetex'] = True


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




def return_or_replicate(data,n):
    if len(data.split(';'))==1:
        return [int(data) for _ in range(n)]
    else:
        return [int(i) for i in data.split(';')]
    
def add_sites(s,array1,array2,array3,array4):
    array1.append(s.data[0][0])
    array3.append(math.sqrt(0.5)*0.529177*(-s.data[1][2]+s.data[1][0]))
    array2.append(s.data[1][1]*0.529177)
    array4.append(math.sqrt(0.5)*0.529177*(-s.data[1][2]-s.data[1][0])/1j)

if __name__=="__main__":

    sns.set_theme("paper","ticks")

    params = {#'text.usetex' : True,
            'font.size' : 11.74983
            #   'font.family' : 'lmodern',
            #   'text.latex.unicode': True,
            }
    plt.rcParams.update(params)
    plt.rcParams["font.family"] = "Nimbus Roman"

    parser = argparse.ArgumentParser(prog="",description="",epilog="")
    parser.add_argument('-f','--files',required=True,help='')
    parser.add_argument('-l','--labels', required=True, help='')
    parser.add_argument('-r','--record_index',required=True,help='')
    parser.add_argument('-s','--site_index',required=True,help='')
    parser.add_argument('-n','--no_of_records',required=True,help='')
    parser.add_argument('-o','--output',required=True,help='')

    args = parser.parse_args()

    args.files=args.files.split(';')
    args.labels=args.labels.split(';')

    if len(args.files)!=len(args.labels):
        raise Exception()
    record_index = return_or_replicate(args.record_index,len(args.files))
    site_index = return_or_replicate(args.site_index,len(args.files))
    no_of_records = return_or_replicate(args.no_of_records,len(args.files))


    order_0 = []
    order_1_0 = []
    order_1_1 = []
    order_1_2 = []
    for i,file_name in enumerate(args.files):
        lines = lib.read_lines(file_name)
        
        DMA_entries = []
        for index in lib.DMA.find(lines,no_of_records[i],file_name):
            DMA_entries.append(lib.DMA.get(lines,index,file_name))

        add_sites(DMA_entries[record_index[i]][site_index[i]],order_0,order_1_0,order_1_1,order_1_2)
        print(order_0[-1])
        print(order_1_0[-1])
        print(order_1_1[-1])
        print(order_1_2[-1])


    
    fig,ax = plt.subplots(figsize=(2.79077651389,0.9*2*0.9),constrained_layout=True)
    l=len(args.labels)
    ax.set_xticks(range(0,l))
    ax.set_xticklabels(args.labels)
    line1, = plt.plot(range(l),order_0,marker='o', linestyle='dashed',label=r"$Q_{00}$",color="red")
    line2, = plt.plot(range(l),order_1_0,marker='o', linestyle='dashed',label=r"$Q_{10}$",color="black")
    line3, = plt.plot(range(l),order_1_1,marker='o', linestyle='dashed',label=r"$Q_{11c}$",color="blue")
    line4, = plt.plot(range(l),order_1_2,marker='o', linestyle='dashed',label=r"$Q_{11s}$",color="green")
    plt.xticks(rotation = 90)
    plt.ylabel(r"$ea_0^l$")
    plt.savefig(args.output+".pdf")
    plt.savefig(args.output+".png",dpi=300)
    plt.clf()
    plt.cla()

    fig,ax = plt.subplots(figsize=(2.79077651389*2,0.25),constrained_layout=True)
    fig.legend((line1,line2,line3,line4), (r"$Q_{00}$", r"$Q_{10}$",r"$Q_{11c}$",r"$Q_{11s}$"),loc="center",ncol=4)
    plt.axis('off')
    # plt.savefig("legend.pdf")


