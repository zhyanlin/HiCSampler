import sys
import click
import glob
import numpy as np
import pandas as pd
from scipy.stats import poisson
import tqdm
from matplotlib import pylab as plt
@click.command()
@click.argument('prefix', type=str,required=True)
@click.argument('obs', type=str,required=True)
@click.argument('output', type=str,required=True)
def pval(prefix,obs,output):
    '''
    prefix: prefix of filenames for all collected RC samples
    obs: observed read counts
    output: pval output
    '''
    obs = pd.read_csv(obs,sep='\t',header=None)
    # print(obs)
    n = np.max(np.max([np.max(obs[1]),np.max(obs[3])])+1)
    obsMat = np.zeros((n,n))
    obsMat[obs[1],obs[3]]=obs[4]
    pMat  = np.zeros((n,n))
    N = 0
    sampleMat = np.zeros((n,n))
    for name in tqdm.tqdm(glob.glob(prefix+"*")):
        rcsample = pd.read_csv(name,sep='\t',header=None)
        rcsample=rcsample[rcsample[1]<n]
        rcsample=rcsample[rcsample[3]<n]
        N = N+1
        sampleMat = sampleMat*0
        sampleMat[rcsample[1],rcsample[3]]=rcsample[4] 
        # pMat += np.multiply((sampleMat>obsMat)*1,(obsMat>0)*1)

        pMat += (sampleMat>=obsMat)*1
    pMat = pMat/N
    pval = pMat[obs[1],obs[3]]
    obs['pval'] = pval
    obs.to_csv(output,sep='\t',header=False, index=False)

    # pp = []
    # for i in range(pMat.shape[0]):
    #     for j in range(i+1,pMat.shape[1]):
    #         pp.append(pMat[i,j])
    # plt.figure()
    # plt.hist(pp,cumulative=True,histtype='step',bins=1000,density=True)
    # plt.plot([0,1],[0,1])
    # plt.show()

    return 0

if __name__ == '__main__':
    pval()