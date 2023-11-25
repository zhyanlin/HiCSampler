# HiCSampler: Posterior inference of Hi-C contact frequency through sampling

Hi-C is one of the most widely used approaches to study three-dimensional genome conformations. Contacts captured by a Hi-C experiment are represented in a contact frequency matrix. Due to the limited sequencing depth and other factors, Hi-C contact frequency matrices are only approximations of the true interaction frequencies, and are further reported without any quantification of uncertainty. Hence downstream analyses based on Hi-C contact maps (e.g TAD and loop annotation) are themselves point estimations.
Here, we present the Hi-C interaction frequency sampler (HiCSampler) that reliably infers the posterior distribution of interaction frequency for a given Hi-C contact map by exploiting dependencies between neighboring loci. Posterior predictive checks demonstrate that HiCSampler is able to infer highly predictive chromosomal interaction frequency. Summary statistics calculated by HiCSampler provide a measurement of the uncertainty for Hi-C experiment, samples inferred by HiCSampler are ready for use by most downstream analyzes tools off the shelf and permit uncertainty measurements in these analyzes without modifications.

## 1. Installation
Please run the following commands to install HiCSampler:
<pre>
git clone https://github.com/zhyanlin/HiCSampler
cd HiCSampler
make
</pre>
You can now run HiCSampler as follows: <code>./HiCSampler </code>

## 2. Example run:

<code>./HiCSampler sampledata/test.RC.txt sampledata/output.txt.gz --bias=./sampledata/test.bias.txt --it=5000 -w=8 --threads=10</code>

This will run HiCsampler for 5000 iterations after the burn-in phase to sample posterior contact maps from read count matrix RC.tsv. It will use a window size of 17x17 to estimate variance in pairwise potentials (17=2*8+1). It outputs summary stats for the posterior distribution to output.txt.gz. In addition, it save samples to folder sampledata.

## 3. Parameters:

--it: number of MCMC interactions.

--threads: number of threads.

--bias: bias vector in ICE normalization.

-w:  window size for variance estimation in pairwise potential. The window size would be 2*w+1 [default 8].

--stepSize: step size for sub-dividing Hi-C contact maps into blocks [default 200].





## 4. Preparing input from data in .[m]cool format
You can use the follow command to convert a HiC data in .[m]cool format into HiCSampler's input format. Both 4.1 and 4.2 output read counts to output.RC.tsv, bias vector to output.bias.tsv

### 4.1 Convert from .mcool file:

<code>
bash ./script/convertfromCool.sh  Your_HiC_File.mcool::/resolutions/resol region resol output 
</code>

### 4.2 Convert from .cool file:

<code>
bash ./script/convertfromCool.sh  Your_HiC_File.cool region resol output 
</code>
