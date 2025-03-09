from fpylll import LLL,BKZ
from g6k import SieverParams
from g6k.siever import Siever
from math import ceil,sqrt,log2
from multiprocessing import Pool
from os import cpu_count, path
import sys
from functools import partial, cache
import math
from g6k import SieverParams


# Some part of the code here was taken from Ludo Pulles in
# https://github.com/ludopulles/AccurateScorePredictionDualSieveAttacks/blob/main/code/utils.py

def change_basis(basis, vector):
	"""
	Puts a vector `vector` specified by coefficients into basis `basis`, i.e.  calculate
	`vector * basis`.
	:param basis: an IntegerMatrix containing the basis for a lattice.
	:param vector: a vector specifying a lattice vector by its coefficients (in terms of `basis`).
	:returns: the same lattice vector as `vector` but now expressed in the canonical basis.
	"""
	#return vector
	return basis.multiply_left(vector)


def norm(vector,q):
	s = 0
	for v in vector:
		vq = v%q
		if(abs(vq) <= abs(vq-q)):
			s += vq**2
		else:
			s+= (vq-q)**2
	return math.sqrt(s)

class Generic_Siever_Param:
	def __init__(self, N, beta_0, beta_1):
		self.N = N
		self.beta_0 = beta_0
		self.beta_1 = beta_1

class Progressive_Siever_Param:
	"""
	def __init__(self,_satrat,_satrad):
		self.satrat = _satrat
		self.satrad = _satrad
		self.l = 0
	"""
	def __init__(self,generic_siever_params,_satrat):
		N = generic_siever_params.N
		beta_0 = generic_siever_params.beta_0
		beta_1 = generic_siever_params.beta_1
		self.satrat = _satrat
		self.satrad = (2*N/_satrat)**(1/beta_1)
		self.l = 0

class Additionnal_Siever_Param:
	def __init__(self,_v = False,_double_db = False,_threads = 1):
		self.v = _v
		self.double_db = _double_db
		self.threads = _threads

class Generic_Siever:
	def __init__(self, generic_siever_params, additional_params):
		self.generic_params = generic_siever_params
		self.add_params = additional_params
	
	def sieve(self,B, bkz=True):
		if(not (self.generic_params.beta_0 == -1)):
			B = BKZ.reduction(B,BKZ.Param(self.generic_params.beta_0))
		#B = BKZ.reduction(B,BKZ.Param(10))
		self.g6k_d = Siever(B, SieverParams(threads=self.add_params.threads, dual_mode=False))
		#B = BKZ.reduction(B,BKZ.Param(self.generic_params.beta_0))
		self.call(self.g6k_d)
		with Pool(self.add_params.threads) as pool:
			self.short_vectors = pool.map(partial(change_basis,B), self.g6k_d.itervalues())

class Progressive_Siever(Generic_Siever):
	#Taken from Accurate Score Prediction
	def __init__(self, generic_siever_params, additional_params, progressive_siever_params):
		Generic_Siever.__init__(self, generic_siever_params, additional_params)
		self.specific_siever_params = progressive_siever_params
	def call(self,g6k):
		self.progressive_sieve(g6k,self.specific_siever_params.l, self.generic_params.beta_1 , saturation_ratio=self.specific_siever_params.satrat, saturation_radius=self.specific_siever_params.satrad, verbose=self.add_params.v,double_db=self.add_params.double_db)
	
	def progressive_sieve(self,g6k, l,  r, saturation_ratio, saturation_radius, verbose, double_db):
		"""
		Sieve in [l, r) progressively. The g6k object will contain a list of short vectors.
		Taking l > 0, will cause G6K to sieve in a projected sublattice of the full basis.

		:param g6k: Siever object used for sieving
		:param l: integer indicating number of basis vectors to skip at the beginning.
		:param r: integer indicating up to where to sieve.
		:param saturation_ratio: ratio to pass on to G6K (i.e. ratio of lattice vectors in the ball of
		radius sqrt(4/3) that should be in the database).
		:param verbose: boolean indicating whether or not to output progress of sieving.
		"""
		verbose = 20
		l = 0
		if verbose:
			print(f"Sieving [{max(l, r-40):3}, {r:3}]", end="", flush=True)
		g6k.initialize_local(l, max(l, r - 40), r)
		g6k(alg="gauss")
		while g6k.l > l:
			# Perform progressive sieving with the `extend_left` operation
			if verbose:
				print(f"\rSieving [{g6k.l:3}, {g6k.r:3}]...", end="", flush=True)
			g6k.extend_left()
			g6k("bgj1" if g6k.r - g6k.l >= 45 else "gauss")
		with g6k.temp_params(saturation_ratio=saturation_ratio, db_size_factor=6):
			g6k(alg="hk3")

		num_dual_vectors = ceil((1 if double_db else 0.5) * saturation_ratio * sqrt(4 / 3)**(r - l))
		if saturation_radius is not None and abs(saturation_radius - sqrt(4/3)) > 1e-6:
			num_dual_vectors = ceil((1 if double_db else 0.5) * saturation_ratio * saturation_radius**(r - l))
			# Note: G6K expects saturation radius passed as a **squared length**, i.e. default is 4/3.
			with g6k.temp_params(saturation_ratio=saturation_ratio, db_size_factor=6,
								db_size_base=saturation_radius, saturation_radius=saturation_radius**2):
				g6k(alg="hk3")

		if verbose:
			print("\rSieving is complete! ", flush=True)
		# Number of dual vectors that is used in a full sieve is (4/3)^{n/2}.
		# Take into account 1/2 since G6K only saves exactly one of the vectors w or -w.
		if len(g6k) < num_dual_vectors:
			print(f"Not enough dual vectors found: {len(g6k)}, expected >= {num_dual_vectors}.")
		assert len(g6k) >= num_dual_vectors
		g6k.resize_db(num_dual_vectors)



