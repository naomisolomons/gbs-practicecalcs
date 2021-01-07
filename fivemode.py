import numpy as np
import sympy as sp
#import scipy as scp
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

    for i in range(len(set)):
        for j in range(len(set)):
            sigma_s[i][j] = covariance[set[i]][set[j]]
            sigma_s[int(len(set))+i][int(len(set))+j] = covariance[int(len(covariance)/2)+set[i]][int(len(covariance)/2)+set[j]]
            sigma_s[i][int(len(set))+j] = covariance[set[i]][int(len(covariance)/2)+set[j]]
            sigma_s[int(len(set))+i][j] = covariance[int(len(covariance)/2)+set[i]][set[j]]

    #Next step is to make further matrices as defined in the notes

    sigma_q = sigma_s + np.identity(2*len(set))

    kappa_s = np.matmul(np.kron(x_matrix, np.identity(len(set))),(np.identity(2*len(set)))-2*np.linalg.inv(sigma_q))

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
            

def Hafnian(kappa_matrix, sigma_q):

    '''This function returns the Hafnian (float) of a given input matrix (square array).'''
    if len(kappa_matrix) == 0:
        hafnian = 1
    else:
        hafnian = 0
    #Probably worth writing in a check that the input is a square matrix and the size is even.
    #First thing we should do is make a list of vectors that represent the different perfectly matching permutations.
    def rec_loop(indices, N, tempmu, allmus, endsize):    
        if len(tempmu) == N:
            #print(tempmu)
            #We have assigned a whole vector. Now we need to add this to our list of vectors, and then reduce our current list
            #appropriately so that we can continue with the algorithm.
            allmus.append(copy.deepcopy(tempmu))
            tempmu.pop(-1)
            tempmu.pop(-1)
            tempmu.pop(-1)
            for i in range(int(N/2)-1):
                if len(allmus) == endsize:
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
                allmus = rec_loop(indices, N, tempmu, allmus, endsize)
        else:
            #This shouldn't happen! But it happened several times when developing this bit of code so I've left this clause here
            #just in case.
            raise ValueError

        return allmus

    def makepair(N):
        '''For the input of a positive even integer N, the output allmus is a list of all the perfectly matching
        pairings of the list of indices 0 to N-1.'''
        #Note that this seems to be by far the slowest bit of the code! It also took the longest to write.
        #It is clearly an inefficient method of finding the PMPs. Needs work!
        if N == 0:
            allmus = []
        elif N == 2:
            allmus = [[0,1]]
        else:
            endsize = sp.factorial2(N - 1)
            listofmodes = list(range(N))
        #We call the recursive function rec_loop.
            allmus = rec_loop(listofmodes, N, [], [], endsize)
        return allmus
    print('before')
    set_M = makepair(len(kappa_matrix))
    print('after')
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

    hafnian = Hafnian(kappa_s_geq_1, sigma_q)
    N = 1
    for p_i in numbers:
        N *= 1/math.factorial(p_i)

    probability = ((2**len(set))*N/np.sqrt(np.linalg.det(sigma_q)))*hafnian
    return probability

#########################################For producing covariance matrices#############################################

def covariance_2mode_squeezed(r, phi):
    c2 = np.cosh(2*r)
    s2 = np.sinh(2*r)
    return [[c2, 0, 0,cmath.rect(s2, phi)],[0,c2,cmath.rect(s2, phi),0],[0,cmath.rect(s2, -phi),c2,0],[cmath.rect(s2, -phi),0,0,c2]]

#set1 = [0,1]
#numbers1 = [1,1]
#covariance1 = covariance_2mode_squeezed(1, math.pi/6)

#print(click_probability(covariance1, displacement0, set1, numbers1))


def random_unitary(n):
    '''This method for generating an nxn Haar random unitary matrix is taken from
    a paper by Francesco Mezzadri, who kindly includes some example code.'''
    z = (np.random.randn(n,n) + 1j*np.random.randn(n,n))/np.sqrt(2.0)
    q,r = np.linalg.qr(z)
    d = np.diagonal(r)
    ph = d/np.abs(d)
    q = np.multiply(q,ph,q)
    return q

def covariance_matrix(unitary1, unitary2, r, mu):
    '''This is the Bloch-Messiah decomp method taken from Oli, Will and Dara's multimode paper.'''
    #That docstring needs work.
    u = unitary1
#    print('u=',u)
    uconj = np.conjugate(u.transpose()) 
    v = unitary2
#    print('v=',v)
    U = np.zeros((2*len(u), 2*len(u)), dtype=np.complex_)
    V = np.zeros((2*len(v), 2*len(v)), dtype=np.complex_)
    vconj = np.conjugate(v.transpose())
    U[:u.shape[0],:u.shape[1]]=u
    U[u.shape[0]:,u.shape[1]:]=uconj
    V[:v.shape[0],:v.shape[1]]=v
    V[v.shape[0]:,v.shape[1]:]=vconj
#    print('U=',U)
#    print('V=',V)
    MD = np.zeros((2*len(r), 2*len(r)), dtype=np.complex_)
    cmat = np.diag(np.cosh(r))
    smat = np.diag(np.sinh(r))
    cmatconj = np.conjugate(cmat.transpose())
    smatconj = np.conjugate(smat.transpose())
    MD[:cmat.shape[0],:cmat.shape[1]]=cmat
    MD[smat.shape[0]:,:smat.shape[1]]=smat
    MD[cmatconj.shape[0]:,cmatconj.shape[1]:]=cmatconj
    MD[:smatconj.shape[0],smatconj.shape[1]:]=smatconj 
#    print('MD=',MD)
    T = np.kron(np.identity(2), np.diag(mu))
#    print('T=',T)
    M_matrix = np.matmul(U, np.matmul(MD, V))
    M_conj = np.conjugate(M_matrix.transpose())
    sigma = np.matmul(M_matrix, np.matmul(T, M_conj))
    print('sigma=',sigma)
    return sigma

#########################################################################################################

numberofmodes = 5
#The way this is set up currently produces a random 5-mode covariance.

unitary1 = random_unitary(numberofmodes)
unitary2 = random_unitary(numberofmodes)
r = [0.5,0.5,0.5,0.5,0.5] #squeezing - would no squeezing be 0?
mu = [1,1,1,1,1] #pure state

random_cov = covariance_matrix(unitary1, unitary2, r, mu)

np.set_printoptions(precision=3)
#print(random_cov - np.conjugate(random_cov.transpose()))

Omega = np.kron(np.identity(numberofmodes), [[0,1],[-1,0]])

list_of_sets = list(product(list(range(2)),repeat=numberofmodes)) #It currently checks for 0,1 photons.
#You could change this by changing range(x) in the line above.
sumofprobs = 0
probabilities = []
for i in range(len(list_of_sets)):
    print(list_of_sets[i])
    prob = click_probability(random_cov, np.zeros(2*numberofmodes), [0,1,2,3,4], list_of_sets[i])
    sumofprobs += prob
    print(prob)
    probabilities.append(prob)
print(sumofprobs)
#print(click_probability(random_cov, np.zeros(10), [0,1,2,3,4], [0,0,0,0,0]))
fig = plt.bar(list(range(len(list_of_sets))), np.real(probabilities))
plt.axhline(linewidth=0.5, color='gray')
plt.xlabel('Click pattern')
plt.ylabel('Probability')
plt.show()


#This is my previous method for generating a random covariance matrix.

#random_A = np.random.rand(5,5) #This is a method from wikipedia for constructing random symplectic matrices
#C_matrix = np.linalg.inv(random_A.transpose())
#D_matrix = np.zeros((10,10))
#D_matrix[:random_A.shape[0],:random_A.shape[1]]=random_A #I took this from stackexchange because there doesn't seem to be a method for direct sum in numpy?
#D_matrix[random_A.shape[0]:,random_A.shape[1]:]=C_matrix #I'm sure I could have worked it out by itself
#random_B = np.random.rand(5,5)
#B_matrix = (random_B + random_B.transpose())/2 #Same here but for symmetrix matrices
#N_matrix = np.identity(10) + np.kron([[0,1],[0,0]], B_matrix)
#Omega_matrix = np.zeros((10,10)) + np.kron([[0,1],[0,0]], np.identity(5)) + np.kron([[0,0],[-1,0]], np.identity(5))
#rand4 = rand.randint(0,3)
#random_symp = np.matmul(np.matmul(np.linalg.matrix_power(Omega_matrix,rand4), D_matrix), N_matrix)
#random_cov = np.matmul(random_symp, random_symp.transpose())
#print(random_cov)


################################TO DO LIST##############################################################

#It would be good to check that everything is input correctly.
#Maybe compare this to the Xanadu method?
