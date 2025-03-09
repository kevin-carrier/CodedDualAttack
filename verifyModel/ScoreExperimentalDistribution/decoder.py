#import sage.all
from interface_polar import *
#from sage.all import Matrix, vector, GF, LinearCode
import galois
import numpy
import copy

class Auxilary_Code:
	def __init__(self,_n,_k,_q):
		self.n = _n
		self.k = _k
		self.q = _q
	def decode(self,x):
		raise NotImplementedError()
	def encode(self,x):
		raise NotImplementedError()
	def unencode(self,x):
		raise NotImplementedError()
	def generator_matrix(self,x):
		raise NotImplementedError()

class Full_Code(Auxilary_Code):
	def __init__(self,_n,_k,_q):
		Auxilary_Code.__init__(self,_n,_k,_q)
	def decode(self,x):
		return [x]
	def dual_syndrome(self,x):
		return x
	def encode(self,x):
		return x
	def unencode(self,x):
		return x
	def generator_matrix(self,x):
		raise NotImplementedError()


class Polar_Code(Auxilary_Code):
	def __init__(self,_n_or_Polar,_k = None,_q = None,_Filter = False,_filter_dlsc = -1):
		if isinstance(_n_or_Polar, Polar_Code):
			o_polar = _n_or_Polar
			_n = o_polar.n
			_k = o_polar.k
			_q = o_polar.q
			Auxilary_Code.__init__(self,_n,_k,_q)
			self.GF = galois.GF(self.q)
			self.Filter = o_polar.Filter
			self.filter_dlsc = o_polar.filter_dlsc
			
			self.C_polar_t = polar_copy(o_polar.C_polar_t)
			self.Mat_Gen = copy.deepcopy(o_polar.Mat_Gen)
			self.Mat_Gen_t = copy.deepcopy(self.Mat_Gen).transpose()
			self.mean_error = o_polar.mean_error 
			self.information = copy.deepcopy(o_polar.information)
		else:
			_n = _n_or_Polar
			Auxilary_Code.__init__(self,_n,_k,_q)
			self.GF = galois.GF(self.q)
			self.Filter = _Filter
			self.filter_dlsc = _filter_dlsc
			print("Generate Polar code")
			while True:
				self.C_polar_t = polar_random(self.q,self.n,self.k)
				#print("Gen Finish")
				nb_fail = 0
				while True:
					#print("Gen Mat")
					_Mat_Gen = []
					for i in range(self.k):
						_Mat_Gen.append(list(polar_random_codeword(self.C_polar_t,self.n)))
					self.Mat_Gen =self.GF(_Mat_Gen)
					#self.code = LinearCode(Mat_Gen)
					if(numpy.linalg.matrix_rank(self.Mat_Gen) == self.k):
						break
					else:
						nb_fail += 1
					if(nb_fail == 5):
						break
				if(nb_fail < 5):
					break
				else:
					polar_free(self.C_polar_t)
			while True:
				#print("INfo")
				I = [int(x) for x in np.random.choice(self.n,self.k,replace=False)]
				I.sort()
				_Mat_Gen_I = [[ M[i] for i in I] for M in _Mat_Gen]
				Mat_Gen_I =self.GF(_Mat_Gen_I)
				if(numpy.linalg.matrix_rank(Mat_Gen_I) == self.k):
					break
			Inv_Mat_Gen_I = numpy.linalg.inv(Mat_Gen_I)
			self.Mat_Gen = numpy.matmul(Inv_Mat_Gen_I, self.Mat_Gen)
			self.mean_error = polar_mean_error(self.C_polar_t)
			if(self.Filter and self.filter_dlsc==-1):
				self.filter_dlsc = self.mean_error
			self.information = I
			self.Mat_Gen_t = copy.deepcopy(self.Mat_Gen).transpose()
	def __del__(self):
		polar_free(self.C_polar_t)
	def decode(self,y):
		if(not self.Filter):
			return [tuple(polar_decode(self.C_polar_t,y))]
		else:
			d = tuple(polar_decode(self.C_polar_t,y))
			n = norm(minus(y,d,self.q),self.q)
			if(n <= self.filter_dlsc):
				return [d]
			else:
				return []
	def generator_matrix(self):
		raise NotImplementedError()
	def encode(self,m):
		raise NotImplementedError()
	def unencode(self,c):
		#print(self.information)
		return tuple([c[i] for i in self.information])
	def check_is_codeword(self,c):
		_m = self.unencode(c)
		m = self.GF([list(_m)])
		cc = [int(numpy.matmul(m,self.Mat_Gen)[0][i]) for i in range(len(c))]
	def dual_syndrome(self,_x):
		x = self.GF([list(_x)])
		return tuple([int(numpy.matmul(x,self.Mat_Gen_t)[0][i]) for i in range(self.k)])
		
