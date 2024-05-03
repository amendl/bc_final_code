import os
from typing_extensions import dataclass_transform
import mbe
import importlib

importlib.reload(mbe)

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

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

def calculate_dma_energy(data1,data2):
    multipole0 = lib.DMA.energy_between_gdma(data1,data2,0)
    multipole1 = lib.DMA.energy_between_gdma(data1,data2,1) + multipole0
    multipole2 = lib.DMA.energy_between_gdma(data1,data2,2) + multipole1
    multipole3 = lib.DMA.energy_between_gdma(data1,data2,3) + multipole2
    multipole4 = lib.DMA.energy_between_gdma(data1,data2,4) + multipole3
    multipole5 = lib.DMA.energy_between_gdma(data1,data2,5) + multipole4
    return [multipole0,multipole1,multipole2,multipole3,multipole4,multipole5]



if __name__=="__main__":

    # load mbe data
    cell=mbe.read_cell('/store/kchfo2/klimj2am/molec_crystals/MeOH/MY_CHECK/input/cell')
    coord, elems=mbe.read_monomers_uni('/store/kchfo2/klimj2am/molec_crystals/MeOH/MY_CHECK/input')
    with open('/home/users/mendla/work/3/data/new_energies_nosym_shell6_ref0_test20tig_AVTZ.dat') as f:
        calc_all=f.readlines()

    # load gdma data
    DMA_entries = [lib.DMA.get_from_gdam_file('/home/users/mendla/work/2/methanol_model_data/MP2_gdma/' + str(file_index)+'/dma4.punch') for file_index in range(4)]      

    # print(DMA_entries[0])  

    # print(coord)
    # print([e.position for e in DMA_entries[0]])
    
    # exit()
    
    calc = calc_all[0:4000]

    print("entries:",len(calc))

    data = np.zeros((len(calc),8))

    dataclass_transform_histogram = np.zeros(int(len(calc)/4))

    for idx,calc_line in enumerate(calc):
        # if idx % 10 != 0: 
        #     continue
        print(idx+1, "out of",len(calc),flush=True)
        # if idx % 4 !=0:
        #     continue

        calc_dir=calc_line.split()[0]
        dir_spl=calc_dir.split('/')
        dir_str=''
        dim_struc=[]
        # XXX
        # print(dir_spl[0],dir_spl[1])

        gdma_data = []
        for dir_act in dir_spl[0:2]:
            vec=mbe.get_dir_vec(dir_act)

            # XXX for checking shift
            # mono=coord[vec[3],:,:].copy()
            # mono_shifted=mbe.mono_shift_vec(mono,vec[0:3],cell)
            # print(mono_shifted)

            # shift monomer accordingly
            gdma_data_part = []
            for e in DMA_entries[vec[3]]:
                new_entry = lib.DMA.Entry(None,None,e)
                new_entry.position = mbe.mono_shift_vec(np.array([new_entry.position]),vec[0:3],cell)[0]/0.529177 
                gdma_data_part.append(new_entry)
            # XXX
            # print([e.position for e in gdma_data_part])

            gdma_data.append(gdma_data_part)
        

        dma_energy = calculate_dma_energy(gdma_data[0],gdma_data[1])
        try: 
            # data_for_histogram[idx] = dma_energy/(float(calc_line.split()[-3])+float(calc_line.split()[-4])) - 1
            data[idx,0]=np.linalg.norm(gdma_data[0][0].position - gdma_data[1][0].position)

            data[idx,1]=float(calc_line.split()[-4])
            data[idx,2]=dma_energy[0]
            data[idx,3]=dma_energy[1]
            data[idx,4]=dma_energy[2]
            data[idx,5]=dma_energy[3]
            data[idx,6]=dma_energy[4]
            data[idx,7]=dma_energy[5]
           
        except:
            pass

        # print(calc_line.split()[-4],dma_energy)
        # print(float(calc_line.split()[-4])/dma_energy)
        # input()
    print("Sorting")
    data = data[data[:,0].argsort()]
    print("Writing")
    with open('test_2.npy', 'wb') as f:
        np.save(f,data)

    print("Done")


    # data = np.histogram(data_for_histogram,np.linspace(_from,_to,_bins+1,endpoint=True))[0]


    # plt.plot(np.linspace(_from,_to - float(_to-_from)/float(_bins),_bins),data, drawstyle='steps-post')
    # # plt.xticks(np.linspace(-0.1,0.1,9,endpoint=True),np.linspace(-0.1,0.1,9,endpoint=True))
    # plt.savefig("mp2_distance.pdf")







    
    # for dir_act in dir_spl[0:2]:
    #     dir_str+=dir_act+'/'
    #     if os.path.isdir(dir_str):
    #         print('dir exists '+dir_str)
    #     else:
    #         print('dir does not exist '+dir_str)
    #         os.mkdir(dir_str)
    #     #get the actual shift from the directory name
    #     vec=mbe.get_dir_vec(dir_act)
    #     #store the monomer structure in the unit cell in the mono object
    #     mono=coord[vec[3],:,:].copy()
    #     #shift the monomer according to the shift
    #     mono_shifted=mbe.mono_shift_vec(mono,vec[0:3],cell)
    #     #create the dimer structure
    #     for i in range(mono.shape[0]):
    #         dim_struc.append([elems[vec[3]][i],' '.join(map(str,mono_shifted[i,:]))])
    # #print the trimer structure to a file
    # with open(dir_str+'dimer.xyz','w+') as f:
    #     f.write(str(len(dim_struc))+'\n')
    #     f.write('generated with py\n')
    #     for lines in dim_struc:
    #         f.write(lines[0]+' '+lines[1]+'\n')




