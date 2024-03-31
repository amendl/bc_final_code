


import random
import math
import os
import numpy as np

import importlib

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

def add_job(alpha,beta,gamma,atoms,geometry,folder,start_script,angles):
    if not os.path.exists(folder):
        os.makedirs(folder)

    rotation_matrix = lib.create_rotation_matrix(alpha,beta,gamma)
    transformed_geometry = np.transpose(rotation_matrix.dot(np.transpose(geometry)))
    with open(folder+"/psi4.com","w") as f:
        f.write("memory 225 gb\n\nmolecule methanol {\nnocom\nnoreorient\n")
        for row_index, row in enumerate(transformed_geometry): 
            f.write(str(atoms[row_index]) + "\t" + str(row[0])+"\t"+str(row[1])+"\t"+str(row[2])+"\n")
        f.write("}\n\nset basis aug-cc-pv5z\nenergy, wfn = energy('scf', return_wfn=True)\nfchk(wfn,'output.fchk')")

    with open(folder+"/data","w") as f:
        f.write("File output.fchk\n\nAngstrom\nMultipoles\n\tswitch 10\n\tLimit 10\n\tPunch dma4.punch\nStart\n\nFinish")
    
    start_script+="cd "+folder+"\n"
    start_script+="sbatch -J "+folder+" ../run.sh\n" 
    start_script+="cd ..\n"

    angles+=str(alpha)+','+str(beta)+','+str(gamma)+"\n"

    return start_script,angles


# Generate a random number between 0 and π


random_number = random.uniform(0, math.pi)

if __name__=="__main__":
    # add_job(0.,0.,0.,geometry.NAMES_1,geometry.GEOMETRY_1,"original","","")
    # exit(0)
    
    start_script = "#! /bin/bash\n"
    angles = ""
    
    
    for i in range(30):
        alpha = random.uniform(0, math.pi)
        beta = random.uniform(0, math.pi)
        gamma = random.uniform(0, math.pi)

        start_script,angles = add_job(alpha,beta,gamma,geometry.NAMES_1,geometry.GEOMETRY_1,str(i),start_script,angles)
    with open('start_all.sh', 'w') as file:
        file.write(start_script)
    with open('angles', 'w') as file:
        file.write(angles)
