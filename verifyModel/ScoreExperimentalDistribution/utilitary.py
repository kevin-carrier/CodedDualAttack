from fpylll import IntegerMatrix
from math import comb
import math

import numpy
import random

from copy import copy

def mod(t,q):
	return tuple([t[i]%q for i in range(len(t))])

def minus(a,b,q):
	return tuple([(a[i] - b[i])%q for i in range(len(a))])

def root_hermite_factor(beta):
	return ((beta/(2*math.pi*math.e))*((math.pi*beta)**(1/beta)))**(1/(2*(beta-1)))

# Volume lattice construction A LWE sample
def volume_lattice(q,m,n):
	return q**(n)
def N_sieve(beta):
	return math.sqrt(4/3)**(beta)

def expected_length_short_vector(vol_lat,d,beta_0,beta_1):
	#print((vol_lat**(1/d))*(N_sieve(beta_1)**(1/beta_1))*(root_hermite_factor(beta_1)**(beta_1-1))*(root_hermite_factor(beta_0)**(d-beta_1)))
	return (vol_lat**(1/d))*(N_sieve(beta_1)**(1/beta_1))*(root_hermite_factor(beta_1)**(beta_1-1))*(root_hermite_factor(beta_0)**(d-beta_1))
	#return (vol_lat**(1/d))*(N_sieve(beta_1)**(1/beta_1))*(root_hermite_factor(beta_1)**(beta_1-1))
	#return 0


def generate_LWE_matrix(m,n,q):
	B = IntegerMatrix.random(n+m, "qary", k=m, q=q)
	A = B.submatrix(0, n, n, m+n)
	return A.transpose()

def construction_A(A,q):
	m = A.nrows
	n = A.ncols
	dim = m+n
	B = IntegerMatrix(dim, dim)
	#At = A.transpose()
	for i in range(m):
		B[i,i] = 1
	for i in range(n):
		B[i+m,i+m] = q
	for i in range(m):
		for j in range(n):
			B[i,j+m] = A[i,j]
	return B



def compute_survival(L, G):
	max_treshold = len(G) - 1
	max_tresh_current = int(max(L)) + 1
	if(max_tresh_current > max_treshold):
		G.extend([0] * (int((max_tresh_current - max_treshold))))
		max_treshold = max_tresh_current
	for l in L:
		if(l >= 0):
			G[int(l)] += 1

def sum_progressive(G):
	for i in range(len(G) - 2, -1, -1):
		G[i] += G[i+1]




class Centered_Binomial:
	def __init__(self,_alpha):
		self.alpha = _alpha
		self.proba = []
		self.value = []
		for i in range(self.alpha+1):
			s = 0
			for j in range(self.alpha+1):
				s += comb(self.alpha,i + j)*comb(self.alpha,j)
			s /= (2**(2*self.alpha))
			self.proba.append(s)
			self.value.append(i)
		self.proba = self.proba + self.proba[1:]
		self.value = self.value + self.value[1:]
		for i in range(self.alpha+1,2*self.alpha + 1):
			self.value[i] = - self.value[i]
	def rand_1(self):
		return int(numpy.random.choice(self.value, 1, p=self.proba)[0])
	def rand(self,n):
		return tuple([self.rand_1() for i in range(n)])


def unif_target(m,q):
	return tuple([random.randrange(0,q) for i in range(m)])
def create_sample(A,alpha_secret, alpha_error,q):
	m = A.nrows
	n = A.ncols
	q = q
	At = copy(A).transpose()
	generator_secret = Centered_Binomial(alpha_secret)
	s = generator_secret.rand(n)
	generator_error = Centered_Binomial(alpha_error)
	e = generator_error.rand(m)
	lc =  At.multiply_left(s)
	b = tuple([lc[i]   + e[i] % q for i in range(m)])
	return b,s,e

def random_matrix(m,n,q):
	return generate_LWE_matrix(m,n,q)

def equalize_and_sum(L_L):
	max_length = max([len(L) for L in L_L])
	for L in L_L:
		L += [0]*(max_length - len(L))
	L_final = [0]*max_length
	for L in L_L:
		for i in range(len(L)):
			L_final[i] += L[i]
	return L_final
def equal(x,y,q):
	for i in range(len(x)):
		if(not x[i]%q == y[i]%q):
			return False
	return True


