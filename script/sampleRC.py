import sys
import click
import numpy as np
import pandas as pd
from scipy.stats import poisson
@click.command()
@click.argument('input', type=str,required=True)
@click.argument('bias', type=str,required=True)
@click.argument('output', type=str,required=True)
def sampleRC(input,bias,output):
    '''
    input: HiCSampler samples in compressed file (.gz)
    out: sampled read counts
    '''
    biasf = pd.read_csv(bias,sep='\t',header=None)
    bias = np.zeros(np.max(biasf[0])+1)
    bias[biasf[0]] = biasf[1]
    df = pd.read_csv(input,sep='\t',header=None)
    df[4]=poisson.rvs(df[4]*bias[1]*bias[3])
    df = df[df[4]>0]
    df.to_csv(output,sep='\t',header=False, index=False)

    return 0

if __name__ == '__main__':
    sampleRC()