import numpy as np
import sympy as sp
from itertools import product
import copy
import math
import cmath

#User inputs

covariance0 = np.zeros((10,10))

displacement0 = np.zeros(10)

#This is a set of indices
set0 = np.zeros((5,), dtype=int)

#This is the set of p_i
numbers0 = [1,1,1,1,1]

########################################

x_matrix = [[0,1],[1,0]]



def kappa_new(covariance, displacement, set, numbers):

    ''' This function gives the kappa matrix with appropriate repetitions of rows and columns according to the 
    p_i values. Inputs:
        covariance: 2N*2N covariance matrix (array)
        displacement: displacement vector
        set: list of modes to be considered
        numbers: list of p_i values
        Outputs:
        kappa_geq_1: the kappa matrix with only the modes we are interested in, repeated an appropriate number of times
        according to the value of the p_i for that mode
        sigma_q: the sigma_q matrix as defined in the notes
    '''

    sigma_s = np.zeros((2*len(set),2*len(set)))

    #This is the covariance matrix but with only the modes we are interested in

    for i in set:
        for j in set:
            sigma_s[i][j] = covariance[i][j]
            sigma_s[int(len(covariance)/2)+i][int(len(covariance)/2)+j] = covariance[int(len(covariance)/2)+i][int(len(covariance)/2)+j]
            sigma_s[i][int(len(covariance)/2)+j] = covariance[i][int(len(covariance)/2)+j]
            sigma_s[int(len(covariance)/2)+i][j] = covariance[int(len(covariance)/2)+i][j]

    #Next step is to make further matrices as defined in the notes

    sigma_q = sigma_s + np.identity(2*len(set))

    kappa_s = np.matmul(np.kron(x_matrix, np.identity(len(set))),(np.identity(2*len(set)))-2*np.linalg.inv(sigma_q))

    #In agreement with doing it by hand up until this point.

    #All this next bit is how we repeat modes as is necessary based on the p_i values

    kappa_s_geq_1 = np.zeros((2*np.sum(numbers),2*np.sum(numbers)))

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

    return kappa_s_geq_1, sigma_q
            

def Hafnian(kappa_matrix):

    '''This function returns the Hafnian (float) of a given input matrix (square array).'''

    hafnian = 0 #placeholder

    #Probably worth writing in a check that the input is a square matrix and the size is even.

    #First thing we should do is make a list of vectors that represent the different perfectly matching permutations.
    #Note that this doesn't seem to work with N = 2. Will need to check this.
    def rec_loop(indices, N, tempmu, allmus):
        if len(tempmu) == N:
            #We have assigned a whole vector. Now we need to add this to our list of vectors, and then reduce our current list
            #appropriately so that we can continue with the algorithm.
            allmus.append(copy.deepcopy(tempmu))
            tempmu.pop(-1)
            tempmu.pop(-1)
            tempmu.pop(-1)
            for i in range(int(N/2)-1):
                if len(allmus) == sp.factorial2(N - 1):
                    break
                elif len(allmus) % sp.factorial2(N - 2*i - 1) == 0:
                    for j in range(N - 2*i -2):
                        tempmu.pop(-1)
                    break
                else:
                    pass

        elif len(tempmu) <N:
            newlist = [index for index in indices if index not in tempmu]
            #The list of indices not yet assigned.
            tempmu.append(newlist[0])
            #Add the first one to our current vector (tempmu).
            for i in range(1,len(newlist)):
                #Now we cycle through the remaining indices. This is how we find all possible pairs.
                tempmu.append(copy.deepcopy(newlist[i]))
                #Once this pair has been assigned, we repeat the method with the remaining indices.
                allmus = rec_loop(indices, N, tempmu, allmus)
        else:
            #This shouldn't happen! But it happened several times when developing this bit of code so I've left this clause here
            #just in case.
            raise ValueError

        return allmus

    def makepair(N):
        '''For the input of a positive even integer N, the output allmus is a list of all the perfectly matching
        pairings of the list of indices 0 to N-1.'''
        listofmodes = list(range(N))
        #We call the recursive function rec_loop.
        allmus = rec_loop(listofmodes, N, [], [])
        return allmus

    set_M = makepair(len(kappa_matrix))
    for mu in set_M:
        pairing = 1
        for k in range(int(len(kappa_matrix)/2)):
            pairing *= kappa_matrix[mu[2*(k+1) -2]][mu[2*(k+1)-1]]
        hafnian += pairing

    return hafnian
        

def click_probability(covariance, displacement, set, numbers):
    '''This is going to give us the probability of seeing a particular pattern of photon numbers (this is the list of integers 'numbers')
    over a particular set of modes (this is the list of integers 'set') with the input covariance matrix (must be size 2N x 2N) and
    displacement vector.'''

    kappa_s_geq_1, sigma_q = kappa_new(covariance, displacement, set, numbers)

    hafnian = Hafnian(kappa_s_geq_1)
    N = 1
    for p_i in numbers:
        N *= 1/math.factorial(p_i)

    probability = ((2**len(set))*N/np.sqrt(np.linalg.det(sigma_q)))*hafnian
    return probability


c2 = 3.7622
s2 = 3.6269

covariance1 = [[c2, 0, 0,s2],[0,c2,s2,0],[0,s2,c2,0],[s2,0,0,c2]]
set1 = [0,1]
numbers1 = [2,2]

print(click_probability(covariance1, displacement0, set1, numbers1))

#Main issues
#2. currently only accepts real numbers - could change it for complex numbers or symbols
#3. I've taken the displacement vector as an input but it doesn't actually use it

#########################################Here's some states to use as examples#############################################