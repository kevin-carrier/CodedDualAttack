import numpy
import math

def dot_product(x,y,q):
	s = 0
	for i in range(len(x)):
		s += x[i]*y[i]
	return s%q

class FFT_sample:
	def __init__(self,_decoded_dual_vectors,_k_fft, _q):
		self.k_fft = _k_fft
		self.q = _q
		self.decoded_dual_vectors = _decoded_dual_vectors
	def init(self,target):
		self.T = numpy.zeros(tuple([self.q for i in range(self.k_fft)]),dtype=numpy.complex128)
		for (decoded, dual_vector) in self.decoded_dual_vectors:
			for i in range(len(decoded)):
				self.T[tuple(decoded)] += math.e**((2j*math.pi/self.q) * dot_product(dual_vector,target,self.q))
	def FFT(self):
		self.T_FFT = numpy.fft.fftn(self.T).real/self.k_fft

class Score_Function:
	def __init__(self,_decoded_dual_vectors,_q):
		self.q = _q
		self.decoded_dual_vectors = _decoded_dual_vectors
	def compute_score(self,error, z):
		self.F = 0
		for (decoded, dual_vector) in self.decoded_dual_vectors:
			self.F += math.cos((2*math.pi/self.q)*(dot_product(dual_vector,error,self.q) - dot_product(decoded,z,self.q)  ))
		return self.F
