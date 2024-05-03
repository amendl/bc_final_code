#! /home/users/mendla/custom_software/env/bin/python

'''

author: amendl@hotmail.com

'''

import matplotlib.pyplot as plt
import numpy as np
import regex
import math


from typing import Final

from scipy.optimize import curve_fit

from sympy import Ynm
from sympy.physics.quantum.cg import Wigner3j
from sympy.functions.combinatorial.factorials import factorial


import quaternionic
import spherical


EPSILON_0 = 8.8541878128E-12
CHARGE    = 1.602176634E-19



def list_table_data(data,summed,start_index,end_index):
    print("start:",start_index," end:",end_index)
    print("========================")
    print("========================")
    print()
    print("start:",data[start_index,0]," end:",data[end_index,0])
    print("n:",end_index-start_index)
    print("SCF:",summed[end_index,1]-summed[start_index,1])
    for i in range(6):
        print(f"Multipoles up to {i}:",summed[end_index,i+2]-summed[start_index,i+2])
    print("========================")
    print("========================")
    
    



def create_rotation_matrix(alpha,beta,gamma):
    a = np.array([[math.cos(alpha),math.sin(alpha),0.],[-math.sin(alpha),math.cos(alpha),0.],[0.,0.,1.]])
    b = np.array([[1.,0.,0.],[0.,math.cos(beta),math.sin(beta)],[0.,-math.sin(beta),math.cos(beta)]])
    c = np.array([[math.cos(gamma),math.sin(gamma),0.],[-math.sin(gamma),math.cos(gamma),0.],[0.,0.,1.]])

    return a.dot(b.dot(c))
    


class PlottingContextManager(object):
    '''
    
    '''
    def __init__(self,path,sufixes=[],show=False):
        self.path       = path
        self.show       = show
        self.sufixes    = sufixes
    def __enter__(self):
        plt.clf()
        fig, ax = plt.subplots(dpi=300)
        return fig,ax
    
    def __exit__(self,exc_type, exc_val, exc_tb):
        plt.legend()
        if not self.sufixes:
            plt.savefig(self.path)
        else:
            for s in self.sufixes:
                plt.savefig(self.path+s)
        if self.show:
            plt.show()

def plot_serie(x,y,label,color=None):
    if color is not None:
        plt.plot(x,y,color=color,linestyle='dashed')
        return plt.scatter(x,y,label=label,color=color)
    else:
        plt.plot(x,y,linestyle='dashed')
        plt.scatter(x,y,label=label)        


def read_lines(path):
    '''

    '''
    error_message = f"molpro output at file \'{path}\' is corrupted."
    f = open(path)
    lines = f.readlines()
    f.close()
    if len(lines)<2:
        raise Exception(error_message)
    if "Molpro calculation terminated" not in lines[len(lines)-1]:
        raise Exception(error_message)
    return lines

def find_in_lines(lines,content):
    '''
    
    '''
    retval = []
    for i,line in enumerate(lines):
        if content in line:
            retval.append(i)
    return retval

def get_value(lines,keyword,characters='='):
    for line in lines:
        if keyword in line:
            return float(line.split(characters)[1].replace(' ','').replace('AU','').replace('mH','').replace("D","E"))
    raise Exception("")

def get_values(lines,keyword,characters='='):
    retval = []
    for line in lines:
        if keyword in line:
            retval.append(float(line.split(characters)[-1].replace(' ','').replace('AU','')))
    return retval


def get_parameters(L,E):
    return E[1]+(E[1]-E[0])/(pow(L/(L-1),3.)-1),(E[1]-E[0])/(pow(1./L,3.)-pow(1./(L-1),3.))

class POP:
    '''
        Class for loading Mulliken population analysis
    '''

    class Error(Exception):
        def __init__(self,message : str):
            super().__init__(message)
    

    KEYWORD         : Final[str] = "PROGRAM * POP (Mulliken population analysis)"
    ERROR_HEADER    : Final[str] = "POP (Mulliken population analysisi) loader"

    @staticmethod
    def find(lines,expected=None,path=None):
        lines = find_in_lines(lines,POP.KEYWORD)
        if expected is not None and len(lines)!=expected:
            raise POP.Error(f"[{POP.ERROR_HEADER}]: Expected {expected} entries but found only {len(lines)} in file \'{path}\'")
        return lines

    @staticmethod
    def get(lines,i,expected_entries,path):
        j=i+1
        retval = []
        while j < len(lines):
            if "Unique atom" in lines[j]:
                for k in range(expected_entries):
                    if len(lines) <= j+k+1 or (not lines[j+k+1]) or lines[j+k+1].isspace():
                        raise POP.Error(f"[{POP.ERROR_HEADER}]: Only {k} instead of {expected_entries} entries loaded from file \'{path}\' at line {j+k+1}")
                    try:
                        retval.append(
                            float(lines[j+k+1].split('+')[-1].split('-')[-1].replace(' ',''))
                            * (-1 if len(lines[j+k+1].split('+'))==1 else 1)
                        )
                    except:
                        raise POP.Error(f"[{POP.ERROR_HEADER}]: Error while reading line {j+k+1} from file \'{path}\'.")
                # print(lines[j+1+expected_entries])
                if lines[j+1+expected_entries]!="" and not lines[j+1+expected_entries].isspace():
                    raise POP.Error(f"[{POP.ERROR_HEADER}]: There might be more entries in file \'{path}\' at line {j+1+expected_entries}")
                return retval
            j+=1
        raise POP.Error(f"[{POP.ERROR_HEADER}]: Cannot find POP data after line {i} in file \'{path}\'.")

    @staticmethod
    def energy(geometry,pop_data):
        energy = 0.
        if len(geometry)!=len(pop_data):
            raise POP.Error(f"[{POP.ERROR_HEADER}]: {len(geometry)} = len(geometry) != len(pop_data) = {len(pop_data)}")
        for j in range(len(geometry)):
            for k in range(len(geometry)):
                if j != k:
                    energy += 0.5*pop_data[j]*pop_data[k]/(np.linalg.norm(geometry[j]-geometry[k])*1e-10)
        return 1./(4*math.pi*EPSILON_0)*energy/(4.3597447222071e-18)*CHARGE*CHARGE

class DMA:
    '''
    
    '''

    class Error(Exception):
        def __init__(self,message : str):
            super().__init__(message)

    class Entry:
        def __init__(self,position,data,orig=None):
            if orig != None:
                self.position = np.copy(orig.position)
                self.l = orig.l
                self.data = []
                for d in orig.data:
                    self.data.append(np.copy(d))
            else:
                if (type(data) != list and type(data) != tuple) or not data:
                    raise Exception("\'data\' should be not empty list or tuple.")
                if type(position)==np.ndarray:
                    if len(position)!=3:
                        raise Exception("\'position\' should be of length 3.")
                    self.position = position
                elif type(position)==list or type(position)==tuple:
                    if len(position)!=3:
                        raise Exception("\'position\' should be of length 3.")
                    self.position = np.array(position)
                else:    
                    raise Exception("\'position\' should be list, tuple or numpy.ndarray.")
                self.l = len(data) -1

                # Add monopole parts
                self.data = [[data[i][0]] for i in range(self.l+1)]
                # Add higher order terms
                for l in range(1,self.l+1):
                    if len(data[l])!=2*l+1:
                        raise Exception(f"Invalid length of data[{i}]. Expected {2*i+1} obtained {len(d)}")
                    for m in range(1,l+1):
                        Qkc = data[l][2*m-1]
                        Qks = data[l][2*m]
                        self.data[l] = [1./math.sqrt(2)*(Qkc-Qks*1.j)] + self.data[l] + [1./math.sqrt(2)*math.pow(-1,m)*(Qkc+Qks*1.j)]
        def rotate(self,q):
            wigner = spherical.Wigner(self.l)
            D           = wigner.D(q)
            new_data = []
            for order,a in enumerate(self.data):
                new_array = np.zeros(2*order + 1,dtype=complex)
                for i1 in range(len(a)):
                    real_index1 = i1 - order
                    for i2 in range(len(a)):
                        real_index2 = i2 - order
                        new_array[i1]+=D[wigner.Dindex(order,real_index1,real_index2)]*a[i2]
                new_data.append(new_array)
            self.data = new_data
            return self
            



            # print(self.data[0])
            # input("")
        
        def __str__(self) -> str:
            retval = f"Multipole at {self.position} of order {self.l}\n"
            for i in range(self.l+1):
                retval += f"Data for order {i}\n"
                retval += f"{self.data[i]}\n"
            return retval
        
        def __repr__(self) -> str:
            return self.__str__()

        @staticmethod
        def energy_old(entry1,entry2):
            '''
                Old implementation based on charge
            '''
            if type(entry1)!=DMA.Entry:
                raise TypeError()
            if type(entry2)!=DMA.Entry:
                raise TypeError()
            return 1./(4*math.pi*EPSILON_0)*entry1.data[0][0]*entry2.data[0][0]*CHARGE*CHARGE/(4.3597447222071e-18)/(np.linalg.norm(entry1.position-entry2.position)*5.29177210903e-11)
        
        @staticmethod
        def energy_of_depth(entry1,entry2,depth):
            '''
                Implementation using multipole expansion (based on spherical harmonics)
                max_depth if 1 works only on charges
            '''
           
            if type(entry1)!=DMA.Entry:
                raise TypeError()
            if type(entry2)!=DMA.Entry:
                raise TypeError()
            if entry1.l < depth:
                raise ValueError()
            if entry2.l < depth:
                raise ValueError()
            energy = 0.
            
            def bracket_factor(l1,l2):
                return math.sqrt(float(factorial(2*l1+2*l2 + 1))/(float(factorial(2*l1))*float(factorial(2*l2))))
            def signum(x):
                if x==0:
                    return 0.
                if x<0:
                    return -1.
                return 1.
                


            # get distance and euler angles for (regular) spherical harmonics
            r_vector = entry1.position - entry2.position
            R = np.linalg.norm(r_vector)
            theta   = math.acos(r_vector[2]/R)
            # if math.isnan(theta):
            #     print("theta")

            # if r_vector[0]*r_vector[0] + r_vector[1]*r_vector[1]==0.:
            #     print("denominator")
            phi     = signum(r_vector[1])*math.acos(r_vector[0]/math.sqrt(r_vector[0]*r_vector[0] + r_vector[1]*r_vector[1]))
            # if math.isnan(phi):
            #     print("phi")
            R_angle     = quaternionic.array.from_spherical_coordinates(theta,phi)
            wigner      = spherical.Wigner(depth)
            Y           = wigner.sYlm(0,R_angle)
            D           = wigner.D(R_angle)


            # sumation over all possible L three indexes
            for l1_index in range(depth+1):
                l2_index = depth - l1_index
                # summation over all m1,m2,m
                for m1_ in range(2*l1_index+1):
                    m1=m1_-l1_index
                    for m2_ in range(2*l2_index+1):
                        m2=m2_-l2_index
                        for m_ in range(2*(l1_index+l2_index)+1):
                            m =m_ - l1_index - l2_index
                            # part = entry1.data[l1_index][l1_index + m1]*entry2.data[l2_index][l2_index + m2]*Y[wigner.Yindex(depth,m)]*spherical.Wigner3j(l1_index,l2_index,depth,m1,m2,m) # XXX
                            # energy +=2*np.real((pow(-1.,l1_index)/pow(2*R,l1_index + l2_index + 1))*bracket_factor(l1_index,l2_index)*part*math.sqrt(4*math.pi/(2*(l1_index+l2_index) + 1)))
                            part = entry1.data[l1_index][l1_index + m1]*entry2.data[l2_index][l2_index + m2]*D[wigner.Dindex(depth,m,0)]
                            energy += np.real((pow(-1.,l2_index)/pow(R,l1_index + l2_index + 1)/pow(0.529177,l1_index+l2_index))*bracket_factor(l1_index,l2_index)*part*spherical.Wigner3j(l1_index,l2_index,depth,m1,m2,m))
            return energy
        @staticmethod
        def energy_of_depth_gdma(entry1,entry2,depth):
            '''
                Implementation using multipole expansion (based on spherical harmonics)
                max_depth if 0 works only on charges
            '''

            if type(entry1)!=DMA.Entry:
                raise TypeError()
            if type(entry2)!=DMA.Entry:
                raise TypeError()
            if entry1.l < depth:
                raise ValueError()
            if entry2.l < depth:
                raise ValueError()
            energy = 0.

            def bracket_factor(l1,l2):
                return math.sqrt(float(factorial(2*l1+2*l2+1))/(float(factorial(2*l1))*float(factorial(2*l2))))
            def signum(x):
                if x==0:
                    return 0.
                if x<0:
                    return -1.
                return 1.



            # get distance and euler angles for (regular) spherical harmonics
            r_vector = entry1.position - entry2.position
            R = np.linalg.norm(r_vector)
            theta   = math.acos(r_vector[2]/R)
            # phi=0
            # if r_vector[0]*r_vector[0] + r_vector[1]*r_vector[1]!=0.:
            #     phi     = signum(r_vector[1])*math.acos(r_vector[0]/math.sqrt(r_vector[0]*r_vector[0] + r_vector[1]*r_vector[1]))
            phi=math.atan2(r_vector[1],r_vector[0])
            # print(theta,phi)

            R_angle     = quaternionic.array.from_spherical_coordinates(theta,phi)
            wigner      = spherical.Wigner(depth)
            Y           = wigner.sYlm(0,R_angle)
            D           = wigner.D(R_angle)


            # sumation over all possible L three indexes
            for l1_index in range(depth+1): 
                l2_index = depth - l1_index
                # summation over all m1,m2,m
                for m1_ in range(2*l1_index+1):
                    m1=m1_-l1_index
                    for m2_ in range(2*l2_index+1):
                        m2=m2_-l2_index
                        for m_ in range(2*(l1_index+l2_index)+1):
                            m =m_ - l1_index - l2_index
                            # part = entry1.data[l1_index][l1_index + m1]*entry2.data[l2_index][l2_index + m2]*Y[wigner.Yindex(depth,m)]*spherical.Wigner3j(l1_index,l2_index,depth,m1,m2,m) # XXX
                            # energy +=2*np.real((pow(-1.,l1_index)/pow(2*R,l1_index + l2_index + 1))*bracket_factor(l1_index,l2_index)*part*math.sqrt(4*math.pi/(2*(l1_index+l2_index) + 1)))
                            part = entry1.data[l1_index][l1_index + m1]*entry2.data[l2_index][l2_index + m2]*D[wigner.Dindex(depth,m,0)]
                            energy += np.real((pow(-1.,l2_index)/pow(R,l1_index + l2_index + 1))*bracket_factor(l1_index,l2_index)*part*spherical.Wigner3j(l1_index,l2_index,depth,m1,m2,m))
            return energy
    

    side_label_regex_gdma         = regex.compile(r'(?<a>[-+]?\d+.\d*) +(?<b>[-+]?\d+.\d*) *+(?<c>[-+]?\d+.\d*) +Rank +(?<d>\d*)')   # use search
    value_regex_gdma               = regex.compile(r'(?<key> *[-+]?\w+.?\w*E?[+-]?\w*)+')
    @staticmethod
    def get_from_gdam_file(path):
        f = open(path)
        lines = f.readlines()
        f.close()

        entries = []


        position = None
        rank     = None
        data_array = None
        data = None


        line_index = 4
        while line_index < len(lines):
            
            # print(lines[line_index].replace('-',' -').replace('+',' +'))
            m = DMA.side_label_regex_gdma.search(lines[line_index].replace('-',' -').replace('+',' +'))
            position = np.array([float(m.group('a')),float(m.group('b')),float(m.group('c'))])
            rank = int(m.group('d'))
            data = []
            for rank_index in range(rank+1):
                data_array = []
                while len(data_array) < 2*rank_index+1:
                    line_index+=1
                    try:
                        m = DMA.value_regex_gdma.match(lines[line_index].replace('+',' +').replace('-',' -'))
                        for v in m.captures('key'):
                            # print(float(v))
                            data_array.append(float(v))
                    except:
                        raise Exception()
                data.append(data_array)
                # print()
            line_index+=2
            # print()
            # print()

            # build dma
            entries.append(DMA.Entry(position=position,data=data,orig=None))
        return entries


        
        # def rotate_expansion(self,phi,theta):
    
    KEYWORD         : Final[str] = "PROGRAM * DMA (Distributed multipole analysis)"
    ERROR_HEADER    : Final[str] = "PROGRAM * DMA (Distributed multipole analysis) loader"

    # .captures('key')
    side_label_regex             = regex.compile(r'(?<a>[-+]?\d+.\d*) +(?<b>[-+]?\d+.\d*) +(?<c>[-+]?\d+.\d*)') # use search
    l_label_regex                = regex.compile(r'^ l( *(?<key>\d+))*') # use match
    value_regex                  = regex.compile(r'Ql\w+[cs]?(?<key> +[-+]?\w+.?\w*E?[+-]?\w*)+') # use match
    
    @staticmethod
    def find(lines,expected=None,path=None):
        lines = find_in_lines(lines,DMA.KEYWORD)
        if expected is not None and len(lines)!=expected:
            raise DMA.Error(f"[{DMA.ERROR_HEADER}]: Expected {expected} entries but found only {len(lines)} in file \'{path}\'")
        return lines
    
    @staticmethod
    def get(lines,i,path):
        j = i + 1
        while j < len(lines):
            if "All multipoles are evaluated up to" in lines[j]:
                try:
                    l = int(lines[j].split('=')[1].replace(' ',''))
                except:
                    raise DMA.Error(f"[{DMA.ERROR_HEADER}]: Error while parsing line {j} in file \'{path}\'.")
                j+=3
                last = False
                retval = []
                while True:
                    # print(lines[j])
                    if not(len(lines)>j+7+2*l):
                        raise DMA.Error(f"[{DMA.ERROR_HEADER}]: At line {j} in \'{path}\' file not long enough")
                    if "Distributed multipoles for" in lines[j]:
                        try:
                            m = DMA.side_label_regex.search(lines[j].replace('-',' -').replace('+',' +'))
                            position = np.array([float(m.group('a')),float(m.group('b')),float(m.group('c'))])
                        except:
                            raise DMA.Error(f"[{DMA.ERROR_HEADER}]: Failed to parse line {j} of \'{path}\'")
                    elif "Total multipoles at the origin" in lines[j]:
                        position = np.array([0.,0.,0.])
                        last = True
                    else:
                        raise DMA.Error(f"[{DMA.ERROR_HEADER}]: Failed to parse line {j} of \'{path}\'")
                    j+=4
                    data = [[] for _ in range(l+1)]
                    for _ in range(1+2*l):
                        try:
                            # print(lines[j])
                            m = DMA.value_regex.match(lines[j])
                            a = m.captures('key')
                            for i,v in enumerate(a):
                                data[l - len(a)+1+i].append(float(v.replace(' ','')))
                        except:
                            raise DMA.Error(f"[{DMA.ERROR_HEADER}]: Failed to parse line {j} of \'{path}\'")
                        j+=1
                    # print(DMA.Entry(position,data).position)
                    # print(DMA.Entry(position,data).data)
                    retval.append(DMA.Entry(position,data))
                    if last:
                        return retval
                    j += 4
                    # input("")
            j+=1
        raise DMA.Error(f"[{DMA.ERROR_HEADER}]: Cannot find DMA data after line {i} in file \'{path}\'.")
    
    @staticmethod
    def energy_old(entries,ignore = [],ignore_single = []):
        retval = 0.
        for j in range(len(entries)):
            for k in range(len(entries)):
                if j==k or (j,k) in ignore or (k,j) in ignore or j in ignore_single or k in ignore_single:
                    pass
                else:
                    retval+=0.5*DMA.Entry.energy_old(entries[j],entries[k])
        return retval
    
    # @staticmethod
    # def energy(entries,ignore = [],ignore_single = [],i=1):
    #     retval = 0.
    #     for j in range(len(entries)):
    #         for k in range(len(entries)):
    #             if j==k or (j,k) in ignore or (k,j) in ignore or j in ignore_single or k in ignore_single:
    #                 pass
    #             else:
    #                 retval+=0.5*DMA.Entry.energy(entries[j],entries[k],i)
    #     return retval
    @staticmethod
    def energy_between(entries1,entries2,depth,ignore1=[],ignore2=[]):
        retval = 0.
        for j in range(len(entries1)):
            for k in range(len(entries2)):
                if j in ignore1 or k in ignore2:
                    pass
                else:
                    retval+=DMA.Entry.energy_of_depth(entries1[j],entries2[k],depth)
        return retval
    @staticmethod
    def energy_between_gdma(entries1,entries2,depth,ignore1=[],ignore2=[]):
        retval = 0.
        for j in range(len(entries1)):
            for k in range(len(entries2)):
                if j in ignore1 or k in ignore2:
                    pass
                else:
                    # print("energy of depth ",depth)
                    retval+=DMA.Entry.energy_of_depth_gdma(entries1[j],entries2[k],depth)
        return retval

    def energy_between_gdma_solving_differences(entries1,entries2,depth,ignore1=[],ignore2=[]):
        retval = 0.
        for j in range(len(entries1)):
            for k in range(len(entries2)):
                if j in ignore1 or k in ignore2:
                    pass
                else:
                    # print("energy of depth ",depth)
                    print("first")
                    print(DMA.Entry.energy_of_depth_gdma(entries2[k],entries1[k],depth))
                    print("second")
                    v = DMA.Entry.energy_of_depth_gdma(entries1[j],entries2[k],depth)
                    retval+=v
                    print(v)
        return retval



    @staticmethod
    def dipole(entries):
        retval = np.zeros((3))
        for i in range(len(entries)-1):
            retval+=entries[i].position*entries[i].data[0][0]
            retval[0]+=0.5*entries[i].data[1][0]
            retval[1]+=0.5*entries[i].data[1][1]
            retval[2]+=0.5*entries[i].data[1][2]
        return retval


class DM:
    '''
    
    '''

    class Error(Exception):
        def __init__(self,message : str):
            super().__init__(message)

    @staticmethod
    def find(lines,path=None,version="1.1",use_debeye=False):
        raise NotImplementedError()
    
    @staticmethod
    def energy(dm1,dm2):
        if type(dm1)!=DM:
            raise DM.Error("")
        if type(dm2)!=DM:
            raise DM.Error("")
        return DMA.Entry.energy(dm1.entry,dm2.entry)

