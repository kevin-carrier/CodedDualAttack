from sage.all import real, sqrt, log, exp, tanh, coth, e, pi, RR, ZZ, sigma, erf, gamma, bessel_J, RealField, binomial, erfinv, floor, ceil, round, find_local_maximum, find_local_minimum, numerical_integral
from estimator.estimator.cost import Cost
from estimator.estimator.lwe_parameters import LWEParameters
from estimator.estimator.lwe import Estimate
from estimator.estimator.reduction import delta as deltaf
from estimator.estimator.reduction import RC, ReductionCost
from estimator.estimator.conf import red_cost_model as red_cost_model_default
from estimator.estimator.util import local_minimum, early_abort_range
from estimator.estimator.io import Logging
from estimator.estimator.nd import NoiseDistribution

from estimator.estimator.schemes import (
    Kyber512,
    Kyber768,
    Kyber1024,
)
import pickle
import multiprocessing
import math

from interface_polar import *

RR = RealField(2048)

img = sqrt(-1)

def frange(a,b,s):
    L=[]
    t = a
    while(t<b):
        L.append(t)
        t += s
    return L

def root_Hermite(b):
    b = RR(b)
    return RR(( (b/(2*pi*e))*((pi*b)**(1/b)) )**(1/(2*(b-1))))

'''
@return the length of the vectors produced by the short vector sampler
@param d: dimension of the lattice
@param V: the volume of the lattice
@param b0: beta_0
@param b1: beta_1
'''
def short_vector_length(d, V, b0, b1):
    d = RR(d)
    V = RR(V)
    b0 = RR(b0)
    b1 = RR(b1)
    return RR( sqrt(4/3) * (V**(1/d)) * (root_Hermite(b1)**(b1-1)) * (root_Hermite(b0)**(d - b1)) )

def Y(alpha, x):
    alpha = RR(alpha)
    x = RR(x)
    if x < RR(10)**(RR(-200)):
        return RR(1.0)
    return RR( gamma(alpha+RR(1)) * bessel_J(alpha, x) / ((x/RR(2))**alpha) )

def get_reduction_cost_model(nn):
    matzov_nns={
        "CN": "list_decoding-naive_classical",
        "CC": "list_decoding-classical",
    }
    if nn in matzov_nns:
        return RC.MATZOV.__class__(nn=matzov_nns[nn])
    elif nn == "C0":
        return RC.ADPS16
    else:
        raise Error("unknown cost model '{}'".format(nn))

#cost_sample(m + nlat, beta0, beta1, N, red_cost_model)
def cost_sample(m, beta0, beta1, N, red_cost_model):
	rho, T, _, beta1 = red_cost_model.short_vectors(
		beta = beta0, N=N, d=m, sieve_dim=beta1
	)
	return rho, T, beta1

def cost_FFT(q,kfft):
	C_mul = RR(1024)
	C_add = RR(160)
	return RR(C_mul*RR(115928) + C_add*RR(240500))*RR(kfft)*RR(q)**(RR(kfft - RR(1)))


def cost_decode(q,nfft,L_size=1):
	C_mul = RR(1024)
	C_add = RR(160)
	L_size = RR(L_size)
	nfft_ = RR(2)**(RR(floor(log(RR(nfft), RR(2)))) + RR(1))
	return RR(3) * L_size * (C_mul * RR(115928) + C_add * RR(240500)) * RR(nfft_) * log(RR(nfft_), RR(2))

def compute_p0(alpha):
	# Centered Binomial distribution
	assert(alpha == 2 or alpha==3)
	prob_B2=[6/16, 4/16, 1/16]
	prob_B3=[20/64, 15/64, 6/64, 1/64]
	if alpha == 2:
		p0 = prob_B2[0]
	else:
		p0 = prob_B3[0]
	return RR(p0)

def survival_normal(T,N):
    res = RR(RR(RR(1) - RR(erf(RR(RR(T)/RR(sqrt(RR(N)))))))/RR(2))
    return RR(log(res, 2))
    

def find_relative_threshold_tmp(alpha, q, m, nenu, nlat, nfft, kfft, beta0, beta1, dlat, dlsc):
	res = RR(beta1 * numerical_integral(lambda t: RR(RR(t**RR(beta1 - 1)) * RR(exp(RR(-alpha*(pi*t*dlat/q)**2)))), 0, 1)[0])
	res *= RR(nfft * numerical_integral(lambda t: RR(RR(t**RR(nfft - 1)) * RR(exp(RR(-alpha*(pi*t*dlsc/q)**2)))), 0, 1)[0])
	return res

def find_relative_threshold(alpha, q, m, nenu, nlat, nfft, kfft, beta0, beta1, dlat, avg_dlsc, sdv_dlsc):
	res = RR(beta1 * numerical_integral(lambda t: RR(RR(t**RR(beta1 - 1)) * RR(exp(RR(-alpha*(pi*t*dlat/q)**2)))), 0, 1)[0])
	res *= RR(exp(RR(-alpha * (pi*avg_dlsc/q)**2 / (1 + 2*alpha*(pi*sdv_dlsc/q)**2))))
	res /= RR(sqrt(1 + 2*alpha*(pi*sdv_dlsc/q)**2))
	return res


def count_s(res):
	i = 0
	d = {}
	for scheme in res:
		d[scheme] = {}
		for nn in res[scheme]:
			d[scheme][nn] = {}
			i += 1
	return i, d



def compute_eta_max(alpha, nenu, nfft):
	res = RR(1.0)
	p0 = RR(compute_p0(alpha))
	
	binom_binom = RR(1.0)
	prob_binom = RR(RR(1.0 - p0)**RR(nenu + nfft))
	for t in range(nenu):
		res -= RR(binom_binom * prob_binom)
		binom_binom *= RR(nenu + nfft - t)
		binom_binom /= RR(t + 1.0)
		prob_binom *= p0
		prob_binom /= RR(1.0 - p0)

	return RR(res)


def compute_eta(R, alpha, nenu, nfft):
	res = RR(1.0)
	p0 = RR(compute_p0(alpha))
	
	binom_binom = RR(1.0)
	prob_binom = RR(RR(1.0 - p0)**RR(nenu + nfft))
	for t in range(nenu):
		res -= RR(binom_binom * prob_binom)
		binom_binom *= RR(nenu + nfft - t)
		binom_binom /= RR(t + 1.0)
		prob_binom *= p0
		prob_binom /= RR(1.0 - p0)

	binom_bis = RR(1.0/binomial(nenu+nfft, nenu))
	for t in range(nenu, nenu+nfft+1):
		#res -= RR( RR( RR( RR(1.0) - RR(binom_bis) )**RR(R) ) * binom_binom * prob_binom)
		res -= RR( RR(exp( RR(R)*RR(log(RR(1.0) - RR(binom_bis))) )) * binom_binom * prob_binom)
		
		binom_bis *= RR(t + 1.0)
		binom_bis /= RR(t + 1.0 - nenu)
		binom_bis = min(RR(1.0), binom_bis)
		binom_binom *= RR(nenu + nfft - t)
		binom_binom /= RR(t + 1.0)
		prob_binom *= p0
		prob_binom /= RR(1.0 - p0)
	return RR(res)

def complexity(alpha, _q, _m, _nenu, _nlat, _nfft, _kfft, red_cost_model, target_proba_false_candidate_global, target_proba_senu, option_dlsc = {'ratio_GV':RR(1.0)}, ret_dict = False):
	d_comp = {}
	q = RR(_q)
	m = RR(_m)
	nenu = RR(_nenu)
	nlat = RR(_nlat)
	nfft = RR(_nfft)
	kfft = RR(_kfft)
	if option_dlsc['experimental_Clsc']:
		assert(option_dlsc['code']['q'] == _q)
		assert(option_dlsc['code']['k'] == _kfft)
		assert(option_dlsc['code']['n'] == _nfft)
		assert(option_dlsc['code']['L_size'] == 1)
		avg_dlsc = option_dlsc['code']['mean_decoding_norm']
		sdv_dlsc = option_dlsc['code']['sigma_decoding_norm']
	else:
		avg_dlsc = RR(option_dlsc['ratio_GV'])*RR( ((RR(q)**(RR(1)-RR(RR(kfft)/RR(nfft)))) * (gamma(nfft/RR(2) + RR(1))**(RR(1)/RR(nfft))) / RR(sqrt(RR(pi)))))
		sdv_dlsc = None
		dlsc = RR(avg_dlsc * (nfft+1)/nfft)
	p0 = RR(compute_p0(alpha))
	R_min = max(RR(1), RR(2 * (p0**(-nenu))))
	R_max = max(R_min, 2^100)
	eta_max = min(compute_eta_max(alpha, _nenu, _nfft),compute_eta(R_max, alpha, _nenu, _nfft))
	if eta_max < target_proba_senu:
		return math.inf, 0, 0

	'''
	R = R_min
	eta = compute_eta(R, alpha, _nenu, _nfft)
	while eta < target_proba_senu:
		R *= 1.1
		eta = compute_eta(R, alpha, _nenu, _nfft)
	'''
	R = RR((R_min + R_max)/2)
	while R_max - R_min > 1:
		eta = compute_eta(R, alpha, _nenu, _nfft)
		if eta < target_proba_senu:
			R_min = R
		else:
			R_max = R
		R = RR((R_min + R_max)/2)
	R = R_max
	eta = compute_eta(R, alpha, _nenu, _nfft)
	
	target_Pwrong = RR(target_proba_false_candidate_global/(R*(q**kfft)))
	
	beta0_inf = 200
	beta0_sup = 1000
	beta0 = (beta0_inf + beta0_sup)//2
	while beta0_sup - beta0_inf > 1:
		rho,_,beta1 = cost_sample(m + nlat, beta0, None, None, red_cost_model)
		avg_dlat = RR(rho * RR(q**RR(nlat/(m+nlat))) * RR(root_Hermite(beta0)**RR(m+nlat-RR(1))))
		#avg_dlat_ = RR(short_vector_length(RR(m+nlat), RR(q)**RR(nlat), RR(beta0), RR(beta1)))
		dlat = RR(avg_dlat * (beta1+RR(1)) / beta1)

		if option_dlsc['experimental_Clsc']:
			relative_threshold = find_relative_threshold(alpha, q, m, nenu, nlat, nfft, kfft, beta0, beta1, dlat, avg_dlsc, sdv_dlsc)
		else:
			relative_threshold = find_relative_threshold_tmp(alpha, q, m, nenu, nlat, nfft, kfft, beta0, beta1, dlat, dlsc)
		
		N = RR(sqrt(4/3))**beta1
		if 2**survival_normal(N*relative_threshold,N) > target_Pwrong :
			beta0_inf = beta0
		else:
			beta0_sup = beta0
		beta0 = (beta0_inf + beta0_sup)//2
	beta0 = beta0_sup
	rho,_,beta1 = cost_sample(m + nlat, beta0, None, None, red_cost_model)
	avg_dlat = RR(rho * RR(q**RR(nlat/(m+nlat))) * RR(root_Hermite(beta0)**RR(m+nlat-RR(1))))
	#avg_dlat_ = RR(short_vector_length(RR(m+nlat), RR(q)**RR(nlat), RR(beta0), RR(beta1)))
	dlat = RR(avg_dlat * (beta1+RR(1)) / beta1)

	if option_dlsc['experimental_Clsc']:
		relative_threshold = find_relative_threshold(alpha, q, m, nenu, nlat, nfft, kfft, beta0, beta1, dlat, avg_dlsc, sdv_dlsc)
	else:
		relative_threshold = find_relative_threshold_tmp(alpha, q, m, nenu, nlat, nfft, kfft, beta0, beta1, dlat, dlsc)
	N = RR(sqrt(4/3))**beta1
	Threshold = RR(N*relative_threshold)
	
	T_FFT = cost_FFT(q,kfft)
	_,T_sample,_ = cost_sample(m + nlat, beta0, beta1, N, red_cost_model)
	
	T_decode = N*cost_decode(q,nfft)
	T_total = T_sample + R*(T_decode + T_FFT)
	if(not ret_dict):
		return float(log(T_total,2)), int(beta0), int(beta1)
	if(ret_dict):
		d_comp['RealField'] = RR
		d_comp['complexity'] = T_total
		d_comp['m'] =  _m
		d_comp['nlat'] = _nlat
		d_comp['nenu'] = _nenu
		d_comp['nfft'] = _nfft
		d_comp['kfft'] = _kfft
		d_comp['beta0'] = beta0
		d_comp['beta1'] = beta1
		d_comp['treshold'] = Threshold
		d_comp['N'] = N
		d_comp['avg_dlsc'] = avg_dlsc
		d_comp['sdv_dlsc'] = sdv_dlsc
		d_comp['avg_dlat'] = avg_dlat
		d_comp['dlat'] = dlat
		d_comp['R'] = R
		d_comp['eta'] = eta
		d_comp['epsilon'] = RR(R*(q**kfft)*(2**survival_normal(Threshold,N)))
		

		d_comp["NT_decode"] = T_decode
		d_comp["T_FFT"] = T_FFT
		d_comp["T_sample"] = T_sample


		d_comp['experimental_Clsc'] = option_dlsc['experimental_Clsc']
		if option_dlsc['experimental_Clsc']:
			d_comp['Clsc_str_repr_code'] = option_dlsc['code']['str_repr_code']
			d_comp['Clsc_nbExperimentDecoding'] =  option_dlsc['code']['aux_nb_decoding_randomWord']
		return d_comp


def gather_statistics_polarCode(q,n,k,L_size, nb_sample_polarization = 50,nb_decodage_meandistance_gen = 50, nb_test_code = 10,nb_decodage_meandistance_stat = 300):
	C_polar = polar_random(q,n,k, nb_test_code = nb_test_code, nb_samples_polarization = nb_sample_polarization, list_size = L_size, nb_test_mean = nb_decodage_meandistance_gen)
	L_d_decode = []
	#d_decode = polar_mean_error(C_polar, list_size = L_size, nb_test = nb_decodage_meandistance)
	stats = polar_stats(C_polar, list_size = L_size, nb_test = nb_decodage_meandistance_stat)
	mean_d_decode = stats[0]
	std_dev_decode = stats[1]
	confidence_decode = stats[2]

	s = polar_str_repr(C_polar)
	polar_free(C_polar)
	res = {}
	res['n'] = n
	res['k'] = k
	res['q'] = q
	res['L_size'] = L_size
	res['mean_decoding_norm'] = mean_d_decode
	res['sigma_decoding_norm'] = std_dev_decode
	res['aux_nb_decoding_randomWord'] = nb_decodage_meandistance_stat

	res['str_repr_code'] = s
	res['aux_confidence_exp'] = confidence_decode
	return res