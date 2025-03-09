from utilitaries import *
"""
with open('optimized_withoutExperimentalPolar_StatisticsPolar.pkl', 'rb') as handle:
    res= pickle.load(handle)

resA = res[0]
resB = res[1]
for scheme in resA:
    for nn in resA[scheme]:
        print(scheme)
        print(nn)
        print(float(resA[scheme][nn]['nfft']))
        print(float(resA[scheme][nn]['kfft']))
        print(float(resB[scheme][nn]['n']))
        print(float(resB[scheme][nn]['k']))
        print(float(resB[scheme][nn]['q']))

        print( "distance")
        print(float(resA[scheme][nn]['mean_dlsc']))
        print(float(resB[scheme][nn]['mean_decoding_norm']))

        #print(float(resA[scheme][nn]['sigma_dlsc']))
        print(float(resB[scheme][nn]['sigma_decoding_norm']))

        print(float(resB[scheme][nn]['aux_nb_decoding_randomWord']))

filename_output3 = 'optimized_withExperimentalPolar.pkl'

"""
from utilitaries import *
with open('optimized_withExperimentalPolar.pkl', 'rb') as handle:
    res= pickle.load(handle)
for scheme in res:
    for nn in res[scheme]:
        param = res[scheme][nn]
        print(scheme)
        print(res)
        print(param)
