from math import comb
import math
from utilitary import *
from sievers import *
import numpy
from copy import copy
#from g6k import BKZ
from fpylll.util import gaussian_heuristic
import FFT_sample
import decoder
import random
from fpylll import FPLLL

from multiprocessing import Process, Pool,Pipe
from functools import partial
from multiprocessing.pool import ThreadPool

def short_vectors(B,beta_0,beta_1,N,verbose = 20):
	# Options for siever
	double_db = False
	sat_ratio = 0.5
	threads = 1
	generic_param = Generic_Siever_Param(N,beta_0,beta_1)
	additionnal_param = Additionnal_Siever_Param(verbose,double_db,threads)

	# Instantiate siever
	progressive_siever_param = Progressive_Siever_Param(generic_param,sat_ratio)
	siever = Progressive_Siever(generic_param,additionnal_param,progressive_siever_param)
	if(verbose != 0):
		print("Begin sieve")
	siever.sieve(B)
	if(verbose != 0):
		print("End sieve")
	return siever.short_vectors

def decode_vectors(C_lsc, A_fft, L_short_vectors, q):
	decoded_vectors = []
	L_norm_lsc_error = []
	m = A_fft.nrows
	L_norm_vector = []
	for short_vector in L_short_vectors:
		L_norm_vector.append(norm(short_vector,q))
		x = short_vector[0:m]
		for maux in C_lsc.decode(mod(A_fft.multiply_left(x),q)):
			#auxilary_code.check_is_codeword(maux)
			L_norm_lsc_error.append(norm(minus(maux,tuple(mod(A_fft.multiply_left(x),q)),q),q))
			decoded_vectors.append((tuple(C_lsc.unencode(maux)),x) )
	return decoded_vectors, L_norm_lsc_error, L_norm_vector 


def algorithm2(C_lsc,A_fft,A_lat,beta_0,beta_1,N,q):
	B_A_lat = construction_A(A_lat,q)
	L_short_vectors = short_vectors(B_A_lat,beta_0,beta_1,N,verbose = 20)
	#C_lsc = decoder.Polar_Code(n_fft,k_fft,q)
	decoded_short_vectors, L_norm_lsc_error, L_norm_vector  = decode_vectors(C_lsc, A_fft, L_short_vectors, q)
	return decoded_short_vectors, L_norm_lsc_error, L_norm_vector 

def score_function_complete(b,decoded_short_vectors,k_fft,q):
	fft = FFT_sample.FFT_sample(decoded_short_vectors,k_fft,q)
	fft.init(b)
	fft.FFT()
	fft_concat = numpy.concatenate(fft.T_FFT,axis=None)
	return fft_concat
	#compute_survival(fft_concat, ret_survival_function)

def survival_uniform_target(A,n_fft, n_lat,k_fft,q,beta_0,beta_1,N,ret_survival_function, C_lsc,number_target):
	m = A.nrows
	n = A.ncols
	assert(n_fft + n_lat == n)
	A_fft = A.submatrix(0, 0, m , 0 + n_fft)
	A_lat = A.submatrix(0,  n_fft, m , n_fft + n_lat)

	#C_lsc = decoder.Polar_Code(n_fft,k_fft,q)

	decoded_short_vectors, L_norm_lsc_error, L_norm_vector = algorithm2(C_lsc,A_fft,A_lat,beta_0,beta_1,N,q)
	for i in range(number_target):
		b = unif_target(m,q)
		fft_concat = score_function_complete(b,decoded_short_vectors,k_fft,q)
		compute_survival(fft_concat, ret_survival_function)
	#return decoded_short_vectors
	return L_norm_lsc_error, L_norm_vector



def _stat_uniform_parallel(m,n,q,n_fft,k_fft,n_lat,beta_0,beta_1,N,option_Clsc, nb_iteration,nb_internal,core,connector):
	FPLLL.set_random_seed(random.randint(0,4294967295))
	ret_survival_function = [0]
	nb_dual = []
	avg_length = 0
	avg_dlsc = 0
	L_norm_dual = []
	L_norm_error = []
	for i in range(nb_iteration):
		A = random_matrix(m,n,q)
		if option_Clsc['change_with_A']:
			C_lsc = decoder.Polar_Code(n_fft,k_fft,q)
		else:
			C_lsc = decoder.Polar_code(option_Clsc['code'])
		L_norm_lsc_error, L_norm_vector = survival_uniform_target(A,n_fft,n_lat,k_fft,q,beta_0,beta_1,N,ret_survival_function,C_lsc,nb_internal)
		nb_dual.append(len(L_norm_vector))
		L_norm_dual += L_norm_vector
		L_norm_error += L_norm_lsc_error
	print(core)
	connector.send([nb_dual, L_norm_dual, L_norm_error, ret_survival_function])

def stat_unif_parallel(m,n,q,n_fft,k_fft,n_lat,beta_0,beta_1,N,option_Clsc,nb_iteration,nb_cores,nb_internal,save_short = False):
	it_per_core = nb_iteration
	fb = partial(_stat_uniform_parallel,m,n,q,n_fft,k_fft,n_lat,beta_0,beta_1,N,option_Clsc,it_per_core,nb_internal)
	total_number_iteration = it_per_core*nb_cores*nb_internal

	L_process = []
	L_pipe = []
	_res = []
	for i in range(nb_cores):
		conn1, conn2 = Pipe()
		p = Process(target=fb, args=(i,conn2,))
		L_process.append(p)
		L_pipe.append((conn1,conn2))
	
	for p in L_process:
		p.start()
	for c1,c2 in L_pipe:
		_res.append(c1.recv())
	for p in L_process:
		p.join()
	L_nb_dual = []
	L_norm_dual = []
	L_norm_error = []
	L_surv = []
	for i in range(nb_cores):
		L_nb_dual += _res[i][0]
		L_norm_dual += _res[i][1]
		L_norm_error += _res[i][2]
		L_surv.append(_res[i][3])
	avg_nb_dual = float(numpy.mean(L_nb_dual))
	ret_survival_function = equalize_and_sum(L_surv)
	sum_progressive(ret_survival_function)


	avg_length = float(numpy.mean(L_norm_dual))
	avg_dlsc = float(numpy.mean(L_norm_error))
	sdv_length = float(numpy.std(L_norm_dual))
	sdv_dlsc = float(numpy.std(L_norm_error))

	diviseur = ((q**k_fft)*total_number_iteration)
	L_ret = []
	for i in range(len(ret_survival_function)):
		ret_survival_function[i] /=diviseur
		L_ret.append(ret_survival_function[i])
	head = "q={:d}\nm={:d}\nn={:d}\nn_fft={:d}\nk_fft={:d}\nn_lat={:d}\nbeta_0={:d}\nbeta_1={:d}\navg_N={:f}\navg_dlat={:f}\nsdv_dlat={:f}\navg_dlsc={:f}\nsdv_dlsc={:f}\nnb_iteration={:d}\nLine i (starting from i=0) correspond to P(F >= i) where F is the score function when the target b is uniform (this is a first approximation of P_wrong(i)). nb_iteration correspond to the number of time we ran the algorithm with each time a different lattice, target".format(q,m,n,n_fft,k_fft,n_lat,beta_0,beta_1,avg_nb_dual,avg_length,sdv_length,avg_dlsc,sdv_dlsc,total_number_iteration)
	name_file = "Pwrong_q{:d}_m{:d}_n{:d}_nfft{:d}_kfft{:d}_nlat{:d}_beta0-{:d}_beta1-{:d}_N{:f}.out".format(q,m,n,n_fft,k_fft,n_lat,beta_0,beta_1,avg_nb_dual)
	numpy.savetxt(name_file, L_ret, delimiter=',',header=head)
	if (save_short):
		numpy.savetxt("NormShortVectors_" + name_file, L_norm_dual, delimiter=',',header=head)
		numpy.savetxt("NormErrorLsc_" + name_file, L_norm_error, delimiter=',',header=head)



# SPECIFIC TARGET
def score_function_target(b,x,decoded_short_vectors,q):
	sf = FFT_sample.Score_Function(decoded_short_vectors,q)
	scoreTarget = sf.compute_score( b, x)
	return scoreTarget

def value_close_target(A,n_fft, n_lat,q,beta_0,beta_1,N,alpha_secret,alpha_error, C_lsc,value_function,number_target): 
	m = A.nrows
	n = A.ncols
	assert(n_fft + n_lat == n)
	A_fft = A.submatrix(0, 0, m , 0 + n_fft)
	A_lat = A.submatrix(0,  n_fft, m , n_fft + n_lat)

	#C_lsc = decoder.Polar_Code(n_fft,k_fft,q)

	decoded_short_vectors, L_norm_lsc_error, L_norm_vector = algorithm2(C_lsc,A_fft,A_lat,beta_0,beta_1,N,q)
	for i in range(number_target):
		_b,_s,_e = create_sample(A,alpha_secret, alpha_error,q)
		b = mod(_b,q)
		s = mod(_s,q)
		e = mod(_e,q)
		s_fft = tuple([s[i]  for i in range(0 ,  n_fft)])
		F_target = score_function_target(b,C_lsc.dual_syndrome(s_fft) ,decoded_short_vectors,q)
		value_function.append(F_target)
	return L_norm_lsc_error, L_norm_vector


def stat_short_vector_sampler(n,m,q,beta0,beta1, N):
	FPLLL.set_random_seed(random.randint(0,4294967295))
	A = random_matrix(m,n,q)
	B_A_lat = construction_A(A,q)
	L_short_vectors = short_vectors(B_A_lat,beta0,beta1,N,verbose = 20)
	L_norm_short_vectors = [norm(x,q) for x in L_short_vectors]

	head = "q={:d}\nn={:d}\nm={:d}\nbeta_0={:d}\nbeta_1={:d}\nN={:f}".format(q,n,m,beta0,beta1,N)
	name_file = "ShortVectors_q={:d}_n={:d}_m={:d}_beta0={:d}_beta1={:d}_N={:d}.out".format(q,n,m,beta0,beta1,N)
	numpy.savetxt(name_file,L_norm_short_vectors, delimiter=',',header=head)

def _stat_target_parallel(m,n,q,n_fft,k_fft,n_lat,beta0,beta1,N,option_Clsc, alpha_secret,alpha_error, nb_iteration,nb_internal,core,connector):
	FPLLL.set_random_seed(random.randint(0,4294967295))
	value_function = []
	nb_dual = []
	avg_length = 0
	avg_dlsc = 0
	L_norm_dual = []
	L_norm_error = []
	for i in range(nb_iteration):
		A = random_matrix(m,n,q)
		if option_Clsc['change_with_A']:
			C_lsc = decoder.Polar_Code(n_fft,k_fft,q)
		else:
			C_lsc = decoder.Polar_code(option_Clsc['code'])
		L_norm_lsc_error, L_norm_vector = value_close_target(A,n_fft, n_lat,q,beta0,beta1,N,alpha_secret,alpha_error, C_lsc,value_function,nb_internal)
		nb_dual.append(len(L_norm_vector))
		L_norm_dual += L_norm_vector
		L_norm_error += L_norm_lsc_error
	print(core)
	connector.send([nb_dual, L_norm_dual, L_norm_error, value_function])

def stat_target_parallel(m,n,q,n_fft,k_fft,n_lat,beta_0,beta_1,N,option_Clsc,alpha_secret,alpha_error,nb_iteration,nb_cores,nb_internal,save_short = False):
	it_per_core = nb_iteration
	fb = partial(_stat_target_parallel,m,n,q,n_fft,k_fft,n_lat,beta_0,beta_1,N,option_Clsc,alpha_secret,alpha_error,it_per_core,nb_internal)
	total_number_iteration = it_per_core*nb_cores*nb_internal
	value_function = []
	L_process = []
	L_pipe = []
	_res = []
	for i in range(nb_cores):
		conn1, conn2 = Pipe()
		p = Process(target=fb, args=(i,conn2,))
		L_process.append(p)
		L_pipe.append((conn1,conn2))
	
	for p in L_process:
		p.start()
	for c1,c2 in L_pipe:
		_res.append(c1.recv())
	for p in L_process:
		p.join()
	L_nb_dual = []
	L_norm_dual = []
	L_norm_error = []
	L_surv = []
	for i in range(nb_cores):
		L_nb_dual += _res[i][0]
		L_norm_dual += _res[i][1]
		L_norm_error += _res[i][2]
		value_function += _res[i][3]
	avg_nb_dual = float(numpy.mean(L_nb_dual))

	avg_length = float(numpy.mean(L_norm_dual))
	avg_dlsc = float(numpy.mean(L_norm_error))
	sdv_length = float(numpy.std(L_norm_dual))
	sdv_dlsc = float(numpy.std(L_norm_error))

	head = "q={:d}\nalpha_secret={:d}\nalpha_error={:d}\nm={:d}\nn={:d}\nn_fft={:d}\nk_fft={:d}\nn_lat={:d}\nbeta_0={:d}\nbeta_1={:d}\navg_N={:f}\navg_dlat={:f}\nsdv_dlat={:f}\navg_dlsc={:f}\nsdv_dlsc={:f}\nnb_iteration={:d}\n Each Line is a value for F(solution)".format(q,alpha_secret,alpha_error,m,n,n_fft,k_fft,n_lat,beta_0,beta_1,avg_nb_dual,avg_length,sdv_length,avg_dlsc,sdv_dlsc,total_number_iteration)
	name_file = "Pgood_q{:d}_alphas{:d}_alphae{:d}_m{:d}_n{:d}_nfft{:d}_kfft{:d}_nlat{:d}_beta0-{:d}_beta1-{:d}_N{:f}.out".format(q,alpha_secret,alpha_error,m,n,n_fft,k_fft,n_lat,beta_0,beta_1,avg_nb_dual)
	numpy.savetxt(name_file, value_function, delimiter=',',header=head)
	if (save_short):
		numpy.savetxt("NormShortVectors_" + name_file, L_norm_dual, delimiter=',',header=head)
		numpy.savetxt("NormErrorLsc_" + name_file, L_norm_error, delimiter=',',header=head)


