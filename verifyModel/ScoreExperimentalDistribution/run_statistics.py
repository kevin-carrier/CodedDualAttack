from Algorithm import *


if __name__ == "__main__":
	for N in [25970]:
		q = 241
		m=40
		n_lat=42
		n_fft = 8
		k_fft = 3
		n = n_lat + n_fft
		beta_0 = 35
		beta_1 = 41

		option_Clsc = {'change_with_A':True}
		nb_cores = 6
		nb_iteration = 10
		nb_internal = 100


		"""
		There will be nb_cores threads running in parallel doing the following:

		do nb_iteration time:
			Choose random A_lat, A_fft
			Choose Clsc according to option_Clsc (If option_Clsc['change_with_A'] = True it chooses a polar code at random. It could also be fixed for all the program at the begining using another option)
			Compute short vectors, decode them
			do nb_internal time:
				choose target (b = uniform or small LWE)
				compute score function
		If b = uniform the probability that F(x) > T is gathered in file Pwrong
		If b = small LWE the value F(solution) is gathered in file Pgood
		"""


		
		# COMPUTE P_WRONG (b = uniform)
		stat_unif_parallel(m,n,q,n_fft,k_fft,n_lat,beta_0,beta_1,N,option_Clsc,nb_iteration,nb_cores,nb_internal,False)
		

		# Compute p_good
		alpha_secret = 2
		alpha_error = 2
		stat_target_parallel(m,n,q,n_fft,k_fft,n_lat,beta_0,beta_1,N,option_Clsc,alpha_secret,alpha_error,nb_iteration,nb_cores,nb_internal,False)
		
		