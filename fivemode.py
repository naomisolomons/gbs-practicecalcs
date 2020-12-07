import numpy as np
import sympy as sp
from itertools import product
import copy
import math
import cmath
import random as rand
#from bristol import ensembles
import matplotlib.pyplot as plt

#User inputs

covariance0 = np.zeros((10,10))

displacement0 = np.zeros(10)
#I've taken the displacement vector as an input but it doesn't actually use it

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

    sigma_s = np.zeros((2*len(set),2*len(set)), dtype=np.complex_)

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

    kappa_s_geq_1 = np.zeros((2*np.sum(numbers),2*np.sum(numbers)), dtype=np.complex_)

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
    print('k mat')
    print(kappa_matrix)

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
    if sum(numbers) < 0:
        print('Need to work on this case.')
        probability = 0
    else:
        kappa_s_geq_1, sigma_q = kappa_new(covariance, displacement, set, numbers)

        hafnian = Hafnian(kappa_s_geq_1)
        N = 1
        for p_i in numbers:
            N *= 1/math.factorial(p_i)

        probability = ((2**len(set))*N/np.sqrt(np.linalg.det(sigma_q)))*hafnian
    return probability

#########################################Here's some examples#############################################

def covariance_2mode_squeezed(r, phi):
    c2 = np.cosh(2*r)
    s2 = np.sinh(2*r)
    return [[c2, 0, 0,cmath.rect(s2, phi)],[0,c2,cmath.rect(s2, phi),0],[0,cmath.rect(s2, -phi),c2,0],[cmath.rect(s2, -phi),0,0,c2]]

set1 = [0,1]
numbers1 = [1,1]
covariance1 = covariance_2mode_squeezed(1, math.pi/6)

#print(click_probability(covariance1, displacement0, set1, numbers1))

#########################################################################################################

#Next step, we are going to try and make some plots with random covariance matrices.

#random_symp = ensembles.Circular().gen_cse(N=5)
#print(random_symp)

random_A = np.random.rand(5,5) #This is a method from wikipedia for constructing random symplectic matrices
C_matrix = np.linalg.inv(random_A.transpose())
D_matrix = np.zeros((10,10))
D_matrix[:random_A.shape[0],:random_A.shape[1]]=random_A #I took this from stackexchange because there doesn't seem to be a method for direct sum in numpy?
D_matrix[random_A.shape[0]:,random_A.shape[1]:]=C_matrix #I'm sure I could have worked it out by itself
random_B = np.random.rand(5,5)
B_matrix = (random_B + random_B.transpose())/2 #Same here but for symmetrix matrices
N_matrix = np.identity(10) + np.kron([[0,1],[0,0]], B_matrix)
Omega_matrix = np.zeros((10,10)) + np.kron([[0,1],[0,0]], np.identity(5)) + np.kron([[0,0],[-1,0]], np.identity(5))
rand4 = rand.randint(0,3)
random_symp = np.matmul(np.matmul(np.linalg.matrix_power(Omega_matrix,rand4), D_matrix), N_matrix)
random_cov = np.matmul(random_symp, random_symp.transpose())
print(random_cov)

list_of_sets = list(product(list(range(3)),repeat=5))
sumofprobs = 0
probabilities = []
#for i in range(len(list_of_sets)):
#    print(list_of_sets[i])
#    sumofprobs += click_probability(random_cov, np.zeros(10), [0,1,2,3,4], list_of_sets[i])
#    print(click_probability(random_cov, np.zeros(10), [0,1,2,3,4], list_of_sets[i]))
#    probabilities.append(click_probability(random_cov, np.zeros(10), [0,1,2,3,4], list_of_sets[i]))
#print(sumofprobs)
print(click_probability(random_cov, np.zeros(10), [0,1,2,3,4], [0,0,0,0,0]))
#plt.bar(list(range(len(list_of_sets))), np.real(probabilities))
#plt.show()

################################TO DO LIST##############################################################

#I need to double check this works when not all modes are selected.
#It would be good to check that everything is input correctly.
#Maybe compare this to the Xanadu method?
