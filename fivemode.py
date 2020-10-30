import numpy as np
import sympy as sp
from itertools import product

#User inputs

covariance0 = np.zeros((10,10))

displacement0 = np.zeros(10)

#This is a set of indices
set0 = np.zeros(5)

#This is the set of p_i
numbers0 = np.zeros(5)

########################################

x_matrix = [[0,1],[1,0]]



def kappa_new(covariance, displacement, set, numbers):

    ''' This function gives the kappa matrix with appropriate repetitions of rows and columns according to the 
    p_i values.
    '''

    sigma_s = np.zeros((2*len(set),2*len(set)))

    for i,j in set:
        sigma_s[i][j] = covariance[i][j]
        sigma_s[2*i][2*j] = covariance[2*i][2*j]
    
    #This is the covariance matrix but with only the modes we are interested in
    
#Is this right??

    sigma_q = sigma_s + np.identity(2*len(set))

    kappa_s = np.matmul(sp.TensorProduct(x_matrix, np.identity(len(set))),(np.identity(2*len(set)))-sigma_q)
    
    kappa_s_geq_1 = np.zeros((2*np.sum(numbers),2*np.sum(numbers)))
    
    #All this next bit is how we repeat modes as is necessary based on the p_i values

    position = 0
    position_temp = 0
    
    numbers2 = list(numbers)
    numbers2.extend(numbers)
    
    for i in range(len(numbers2)):
        for j in range(len(numbers2)):
            for k,l in product(range(numbers2[i]),range(numbers2[j])):
                kappa_s_geq_1[position+k][position_temp+l] = kappa_s[i][j]
            position_temp += numbers2[j]
        position_temp = 0
        position += numbers2[i]

    return kappa_s_geq_1
            

def Hafnian(kappa_matrix):

    hafnian = 0 #placeholder

    #The set size is the number of modes we are interested in. If this is S, the size of kappa_matrix should be 2S x 2S
    set_size = len(kappa_matrix)

    #Probably worth writing in a check that the input is a square matrix and the size is even.

    #First thing we should do is make a list of vectors that represent the different perfectly matching permutations.

    indices = list(range(set_size))
    def make_pair(indices):
        
        all_mus = []
        mu = []
        if len(indices)>0:
            mu.append(indices[0])
            indices.pop(0)
            for i in range(len(indices)):
                mu.append(indices[i])
                indices.pop(i)
                make_pair(indices, full_list)
        else:
            all_mus.append(mu)
        
        return all_mus




    return hafnian
        

def make_pair(set_size):
    all_mus = []
    temp_mu = []
    indices = list(range(set_size))
    
    def loop_rec(all_mus, temp_mu, indices):
        if len(indices)>0:
            for i in range(len(indices)):
                temp_mu.append(indices[0])
                indices.pop(0)
                print(temp_mu)
                temp_mu.append(indices[i])
                indices.pop(i)
                loop_rec(all_mus, temp_mu, indices)
        else:
            all_mus.append(temp_mu)
    
    all_mus = loop_rec(all_mus, temp_mu, indices)
    return all_mus

#indices = list(range(6))

check = make_pair(6)
print(check)

########################################################Old code used for testing#################################################
#numbers = [3,1]

#numbers2 = list(numbers)
#numbers2.extend(numbers)

#kappa_s_geq_1 = np.zeros((2*np.sum(numbers),2*np.sum(numbers)))

#kappa_s = [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]

#position = 0
#position_temp = 0
    
#for i in range(len(numbers2)):
 #   for j in range(len(numbers2)):
       # position_temp = position
 #       for k,l in product(range(numbers2[i]),range(numbers2[j])):
  #          kappa_s_geq_1[position+k][position_temp+l] = kappa_s[i][j]
            #kappa_s_geq_1[sum(numbers)+position+i][sum(numbers)+position_temp+j] = kappa_s[len(numbers)+i][len(numbers)+j]
            #print(kappa_s_geq_1)
#        position_temp += numbers2[j]
 #   position_temp = 0
 #   position += numbers2[i]

#print(kappa_s)
#print(kappa_s_geq_1)