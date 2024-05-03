

import numpy as np

NAMES_1     = ['C','H','O','H','H','H']
GEOMETRY_1  = np.array(
    [
        [0.90001,2.35489,0.93548],
        [0.14082,2.42126,1.66879],
        [2.16840,2.20452,1.52603],
        [2.40716,3.03203,2.05362],
        [0.73287,1.45081,0.35202],
        [0.89513,3.16894,0.24385]
    ]
)
NAMES_2     = ['C','H','O','H','H','H']
GEOMETRY_2  = np.array(
    [
        [1.53639, 2.28621, 5.36903],
        [2.29558, 2.21984, 6.10234],
        [0.268, 2.43658, 5.95958],
        [0.02924, 1.60907, 6.48717],
        [1.70353, 3.19029, 4.78557],
        [1.54127, 1.47216, 4.67739]
    ]
)
NAMES_3     = ['C','H','H','O','H','H']
GEOMETRY_3  = np.array(
    [
        [3.97279, 0.03434, 3.49807],
        [3.97767, 0.84839, 4.18970],
        [4.13993,-0.86974, 4.08153],
        [2.70440,-0.11603, 2.90752],
        [2.46564, 0.71148, 2.37993],
        [4.73198, 0.10071, 2.76476]
    ]
)


def cross_product_matrix(vector):
    return np.array(
        [
            [0.,-vector[2],vector[1]],
            [vector[2],0.,-vector[0]],
            [-vector[1],vector[0],0.]
        ]
    )

def create_rotation_matrix(line,angle):
    '''

    '''
    line = line / np.linalg.norm(line)
    return np.cos(angle)*np.identity(3) + np.sin(angle)*cross_product_matrix(line) +(1-np.cos(angle))*np.outer(line,line)

def rotate_molecule(positions,rotation_matrix,center_index):
    '''
    
    '''
    return np.transpose(np.matmul(rotation_matrix,np.transpose(positions-positions[center_index])))+positions[center_index]
def move_molecule(positions,vector):
    '''
    
    '''
    return positions+vector

class LineScanOHBond01:
    '''
        scans 10 positions starting at original position and ending at 4 times the original distance of monomers
    '''
    @staticmethod
    def generator():
        vector = GEOMETRY_3[4]-GEOMETRY_1[2]
        n,division = 10,3
        for i in range(n):
            yield i
    @staticmethod
    def transform(id):
        vector = GEOMETRY_3[4]-GEOMETRY_1[2]
        n,division = 10,3
        return (vector*float(id)/float(division),np.identity(3))
    @staticmethod
    def plot_position(id):
        n,division = 10,3
        return (float(id)/float(division) + 1.)*np.linalg.norm(GEOMETRY_3[4]-GEOMETRY_1[2])
    @staticmethod
    def file_index(id):
        return id
    
class LineScanHF01:
    '''
        scans 10 positions starting at original position and ending at 4 times the original distance of monomers
    '''
    @staticmethod
    def generator():
        n,division = 10,3
        for i in range(n):
            yield i
    @staticmethod
    def transform(id):
        vector = np.array([0.917*2,0.,0.])
        n,division = 10,3
        return (vector*float(id)/float(division),np.identity(3))
    @staticmethod
    def plot_position(id):
        n,division = 10,3
        return (float(id)/float(division) + 1.)*np.linalg.norm(GEOMETRY_3[4]-GEOMETRY_1[2])
    @staticmethod
    def file_index(id):
        return id


class LineScanCO201:
    '''
        scans 10 positions starting at original position and ending at 4 times the original distance of monomers
    '''
    @staticmethod
    def generator():
        n,division = 10,3
        for i in range(n):
            yield i
    @staticmethod
    def transform(id):
        vector = np.array([2.,0.,0.])
        n,division = 10,3
        return (vector*float(id)/float(division),np.identity(3))
    @staticmethod
    def plot_position(id):
        n,division = 10,3
        return (float(id)/float(division) + 1.)*np.linalg.norm(np.array([2.,0.,0.]))
    @staticmethod
    def file_index(id):
        return id

class LineScanOHBond01Close:
    @staticmethod
    def generator():
        for i in range(5):
            yield i
    @staticmethod
    def transform(id):
        vector = GEOMETRY_3[4]-GEOMETRY_1[2]
        return(-vector*0.5 + vector*float(id)/5.,np.identity(3))
    @staticmethod
    def plot_position(id):
        return 0.5+float(id)/5.
    @staticmethod
    def file_index(id):
        return id

    
class LineScanOHBond02:
    '''
        scans 10 positions starting at original position and ending at 4 times the original distance of monomers with rotations
    '''
    @staticmethod
    def generator():
        n,division = 10,3
        for i in range(n):
            for j in range(3):
                yield i,j
    @staticmethod
    def transform(id):
        vector = GEOMETRY_3[4]-GEOMETRY_1[2]
        n,division = 10,3
        return (vector*float(id[0])/float(division),create_rotation_matrix(vector,2*id[1]*3.141592/3.))
    @staticmethod
    def plot_position(id):
        n,division = 10,3
        return float(id[0])/float(division) + 1.
    @staticmethod
    def file_index(id):
        return id[0]*3 + id[1]
    

if __name__=="__main__":
    print(np.linalg.norm(GEOMETRY_3[4]-GEOMETRY_1[2]))
