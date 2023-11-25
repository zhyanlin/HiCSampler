import sys
import click
import glob
import numpy as np
import pandas as pd
from scipy.stats import poisson
import tqdm
from matplotlib import pylab as plt
@click.command()
@click.argument('posteriorIF', type=str,required=True)
@click.argument('bias', type=str,required=True)
@click.argument('obs', type=str,required=True)
def evidence(posteriorif,bias,obs):
    '''
    posteriorIF: prefix of filenames for all collected RC samples
    obs: observed read counts
    bias: bias file
    '''
    obs = pd.read_csv(obs,sep='\t',header=None)
    biasf = pd.read_csv(bias,sep='\t',header=None)
    bias = np.zeros(np.max(biasf[0])+1)
    bias[biasf[0]] = biasf[1]
    bias[bias<0.1] = 1
    IF = pd.read_csv(posteriorif,sep='\t',header=None)
    IF[4]=IF[4]*bias[1]*bias[3]

    n = np.max(np.max([np.max(IF[1]),np.max(IF[3])])+1)
    obsMat = np.zeros((n,n))
    obsMat[obs[1],obs[3]]=obs[4]
    IFMat= np.zeros((n,n))
    IFMat[IF[1],IF[3]]=IF[4]


    # ll=0
    # for i in tqdm.tqdm(range(obsMat.shape[0]-10)):
    #     for j in range(i+1,obsMat.shape[1]-10):
    #         # if obsMat[i,j]>-1:
    #         ijll=np.log(poisson.pmf(obsMat[i,j], np.clip(IFMat[i,j],a_min=1e-5,a_max=1e10)))
    #         # print('\t'.join([str(x) for x in [i,j,ijll]]))
    #         ll+=ijll
    # print("ll=\t"+str(ll))
    llmaps=np.log(poisson.pmf(obsMat, np.clip(IFMat,a_min=1e-5,a_max=1e10)))
    ll=llmaps[0:obsMat.shape[0]-10,0:obsMat.shape[0]-10][np.triu_indices(obsMat.shape[0]-10, k=1)].sum()
    print("ll=\t"+str(ll))
    return 0

if __name__ == '__main__':
    evidence()