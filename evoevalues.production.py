import sys
import numpy as np
import pandas as pd
import joblib
from statsmodels.distributions import ECDF



evos=[pd.read_csv(sys.argv[1]+'.evo'+str(i),sep='\t',comment='N',header=None,index_col=[0,1],low_memory=False) for i in range(1,4)]
evos=pd.concat(evos,axis=1)

print evos.shape

def impute(x):
    defaultvals=[0.115,-0.033,1.909,-0.200]
    for i in range(x.shape[1]):
        x[np.isnan(x[:,i]),i]=defaultvals[i]
    return x

evos.columns=np.arange(4)
evos.iloc[:,:]=impute(evos.values)

evos.to_csv(sys.argv[1]+'.evoall',header=False,float_format ='%.6e')

allevos_ecdf=joblib.load('./resources/ecdfs/allevos_ecdf.pkl')

evos_e=np.zeros((evos.shape[0],4))
evos_e[:,0]=1-allevos_ecdf[0](np.asarray(evos.iloc[:,0]))
evos_e[:,1]=1-allevos_ecdf[3](np.asarray(evos.iloc[:,1]))
evos_e[:,2]=1-allevos_ecdf[6](np.asarray(evos.iloc[:,2]))
evos_e[:,3]=1-allevos_ecdf[7](np.asarray(evos.iloc[:,3]))
evos_e[evos_e==0]=1e-6

evos.iloc[:,:]=evos_e

evos.to_csv(sys.argv[1]+'.evo.evalues',header=False,float_format ='%.6e')

