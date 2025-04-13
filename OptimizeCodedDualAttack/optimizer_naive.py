from utilitaries import *
from functools import partial
import math

#Optimizer

def optimize_init(alpha, q, n, red_cost_model, target_proba_false_candidate_global, target_proba_senu, option_dlsc):
	L_param_best = [0,0,0,0,0,0]
	min_Comp = math.inf
	m = n
	for m in range(n-100, n+1, 20):
		for nenu in range(5, 45, 10):
			for nfft in range(10, 101, 10):
				nlat = n - nfft - nenu
				for kfft in range(1,nfft//2+2,5):
					Comp, beta1, beta0 = complexity(alpha, q, m, nenu, nlat, nfft, kfft, red_cost_model, target_proba_false_candidate_global, target_proba_senu, option_dlsc = option_dlsc)
					if (Comp < min_Comp):
						min_Comp = Comp
						L_param_best = [m,nenu,nfft,kfft,beta0,beta1]
						print("BETTER ------------------------------- ")
						print(float(min_Comp))
						print(L_param_best)
		
	return L_param_best, min_Comp

def optimize_from(alpha, q, n, red_cost_model, target_proba_false_candidate_global, target_proba_senu, L_param, option_dlsc, lock_nfft_kfft):
	if(not lock_nfft_kfft):
		around_m = 1
		around_nenu = 1
		around_nfft = 1
		around_kfft = 1
	else:
		around_m = 3
		around_nenu = 3
		around_nfft = 0
		around_kfft = 0
							 
	L_param_best = L_param.copy()
	[m_best,nenu_best,nfft_best,kfft_best,beta0_best,beta1_best] = L_param
	nlat_best = n - nenu_best - nfft_best
	min_Comp, beta0, beta1 = complexity(alpha, q, m_best, nenu_best, nlat_best, nfft_best, kfft_best, red_cost_model, target_proba_false_candidate_global, target_proba_senu, option_dlsc = option_dlsc)
	for m in range(max(m_best - around_m,1) , min(m_best + around_m,n) + 1):
		for nenu in range(max(nenu_best - around_nenu,1),min(nenu_best + around_nenu,100)  +1):
			for nfft in range(max(nfft_best - around_nfft,1),min(nfft_best + around_nfft,200)+1):
				nlat = n - nfft - nenu
				for kfft in range(max(kfft_best - around_kfft,1),min(kfft_best + around_kfft,nfft - 1)+1):
					Comp, beta0, beta1 = complexity(alpha, q, m, nenu, nlat, nfft, kfft, red_cost_model, target_proba_false_candidate_global, target_proba_senu, option_dlsc = option_dlsc)
					if (Comp < min_Comp):
						min_Comp = Comp
						L_param_best = [m,nenu,nfft,kfft,beta0,beta1]
						print("BETTER ------------------------------- ")
						print(float(min_Comp))
						print(L_param_best)

	return L_param_best, min_Comp

def optimize_1(alpha, q, n, red_cost_model, target_proba_false_candidate_global, target_proba_senu, L_param, nb_It, option_dlsc, lock_nfft_kfft): 
	#if not lock_nfft_kfft:
	#	L_param, Comp = optimize_init(alpha, q, n, red_cost_model, target_proba_false_candidate_global, target_proba_senu, option_dlsc)

	[m_best,nenu_best,nfft_best,kfft_best,beta0_best,beta1_best] = L_param
	nlat_best = n - nenu_best - nfft_best
	Comp, beta0_best, beta1_best = complexity(alpha, q, m_best, nenu_best, nlat_best, nfft_best, kfft_best, red_cost_model, target_proba_false_candidate_global, target_proba_senu, option_dlsc)
	sComp_base = Comp
	
	for i in range(nb_It):
		save_comp = Comp
		L_param, Comp = optimize_from(alpha, q, n, red_cost_model, target_proba_false_candidate_global, target_proba_senu, L_param, option_dlsc, lock_nfft_kfft)
		if(Comp == save_comp):
			if(Comp == sComp_base):
				print("END NOT BETTER ------------------------------ " + str(n))
			return L_param, Comp
	return L_param, Comp


def optimize_(res_starting_point, target_proba_false_candidate_global, target_proba_senu, nb_iteration_optimiser, option_dlsc,ii):
	i = 0
	for scheme in res_starting_point:
		for nn in res_starting_point[scheme]:
			if(i != ii):
				i += 1
			else:
				parameters = res_starting_point[scheme][nn]
				
				red_cost_model = get_reduction_cost_model(nn)
				
				cst = scheme.normalize()
				m = parameters['m']
				nenu = parameters['nenu'] 
				nfft = parameters['nfft'] 
				kfft =parameters['kfft'] 
				beta0 = parameters['beta0'] 
				beta1 =parameters['beta1']
				L_size = 1
				n = cst.n
				q=3329
				nlat = n - nenu - nfft
				if(n == 512):
					alpha = 3
				else:
					alpha = 2
				L_param = [m,nenu,nfft,kfft,beta0,beta1]
				
				nb_It = nb_iteration_optimiser

				option_dlsc_ = option_dlsc[scheme][nn]
				if option_dlsc_['experimental_Clsc']:
					lock_nfft_kfft = True
				else:
					lock_nfft_kfft = False
				L_best,Comp_best = optimize_1(alpha, q, n, red_cost_model, target_proba_false_candidate_global, target_proba_senu, L_param, nb_It, option_dlsc = option_dlsc_, lock_nfft_kfft = lock_nfft_kfft)
				print(scheme)
				print(nn)
				print(L_best)
				print(Comp_best)
				return [scheme, nn, L_best, Comp_best]

def optimize(res_starting_point, nb_iteration_optimiser, nb_core, option_dlsc):
	target_proba_false_candidate_global = RR(0.05)
	target_proba_senu = RR(0.6)
	q = 3329
	#option_dlsc = {'ratio_GV': 1.1}
	#short_description += ", with " + str(nb_iteration_optimiser) + " optimizer iteration"
	nbS, new_comp = count_s(res_starting_point)
	func_ = partial(optimize_, res_starting_point,target_proba_false_candidate_global,target_proba_senu, nb_iteration_optimiser,option_dlsc)
	pool = multiprocessing.Pool(processes = nb_core)
	LL = pool.map(func_, range(nbS))
	pool.close()
	pool.join()
	for [scheme, nn, L_best, Comp_best] in LL:
		[m,nenu,nfft,kfft,beta0,beta1] = L_best
		red_cost_model = get_reduction_cost_model(nn)
		cst = scheme.normalize()
		n = cst.n
		nlat = n - nenu - nfft
		if(n == 512):
			alpha = 3
		else:
			alpha = 2
		option_dlsc_ = option_dlsc[scheme][nn]
		d_comp = complexity(alpha, q, m, nenu, nlat, nfft, kfft, red_cost_model, target_proba_false_candidate_global, target_proba_senu, option_dlsc_, ret_dict = True)
		new_comp[scheme][nn] = d_comp
	return new_comp




def optimize_from_starting_parameter_without_experimental_polar_code(filename_starting_point,filename_output):
	with open(filename_starting_point, 'rb') as handle:
		res_starting_point= pickle.load(handle)

	nb_iteration_optimiser = 700

	option_dlsc = {}
	for scheme in res_starting_point:
		option_dlsc[scheme] = {}
		for nn in res_starting_point[scheme]:
			option_dlsc[scheme][nn] = {'experimental_Clsc':False, 'ratio_GV':RR(1.0)}

	res = optimize(res_starting_point, nb_iteration_optimiser, nb_core = 7, option_dlsc = option_dlsc)
	with open(filename_output, 'wb') as f:
		pickle.dump(res, f)


def generate_data_polar_code_(res, ii):
	i = 0
	for scheme in res:
		for nn in res[scheme]:
			if(i != ii):
				i += 1
			else:
				L_size = 1
				q = 3329
				n = res[scheme][nn]['nfft']
				k = res[scheme][nn]['kfft']

				nb_sample_polarization = 100
				nb_decodage_meandistance_gen = 100
				nb_test_code = 10
				nb_decodage_meandistance_stat = 1000
				
				d = gather_statistics_polarCode(q,n,k,L_size, nb_sample_polarization = nb_sample_polarization,nb_decodage_meandistance_gen = nb_decodage_meandistance_gen, nb_test_code = nb_test_code,nb_decodage_meandistance_stat = nb_decodage_meandistance_stat)
				print(d)
				return [scheme, nn , d]


def generate_data_polar_code(filename, output):
	with open(filename, 'rb') as handle:
		res = pickle.load(handle)
	nbS, resB = count_s(res)
	func_ = partial(generate_data_polar_code_, res)
	pool = multiprocessing.Pool(processes = 7)
	LL = pool.map(func_, range(nbS))
	pool.close()
	pool.join()
	for [scheme, nn , d] in LL:
		resB[scheme][nn] = d
	resF = {}
	resF[0] = res
	resF[1] = resB
	with open(output, 'wb') as f:
		pickle.dump(resF, f)

def optimize_from_starting_parameter_with_experimental_polar_code(filename_starting_point,filename_output):
	with open(filename_starting_point, 'rb') as handle:
		res_starting_point= pickle.load(handle)
	
	res = res_starting_point[0]
	Clsc = res_starting_point[1]


	nb_iteration_optimiser = 700

	option_dlsc = {}
	for scheme in res:
		option_dlsc[scheme] = {}
		for nn in res[scheme]:
			cc = Clsc[scheme][nn].copy()
			option_dlsc[scheme][nn] = {'experimental_Clsc':True, 'code':cc}

	res = optimize(res, nb_iteration_optimiser, nb_core = 7, option_dlsc = option_dlsc)

	with open(filename_output, 'wb') as f:
		pickle.dump(res, f)

if __name__ == "__main__": 
	# Supposed to run in a night maximum
	filename_starting_point = 'start_parameter.pkl'
	filename_output1 = 'optimized_withoutExperimentalPolar.pkl'
	filename_output2 = 'optimized_withoutExperimentalPolar_StatisticsPolar.pkl'
	filename_output3 = 'optimized_withExperimentalPolar.pkl'


	optimize_from_starting_parameter_without_experimental_polar_code(filename_starting_point,filename_output1)
	generate_data_polar_code(filename_output1,filename_output2)
	optimize_from_starting_parameter_with_experimental_polar_code(filename_output2, filename_output3)