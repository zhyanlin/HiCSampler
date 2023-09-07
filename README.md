# HiCSampler: Posterior inference of Hi-C contact frequency through sampling

### Parameters:
--it: number of MCMC interactions.

--threads: number of threads.

--bias: bias vector in ICE normalization (https://github.com/hiclib/iced).

--stepSize: step size for sub-dividing Hi-C contact maps into blocks [default 200].



### Example run:

<code>HiCSampler sampledata/test.RC.txt sampledata/output.txt.gz --bias=sampledata/bias.txt --it=5000 --threads=10</code>

This will run HiCsampler for 5000 iterations after the burn-in phase to sample posterior contact maps from read count matrix RC.tsv. It outputs summary stats for the posterior distribution to output.txt.gz. In addition, it save samples to folder sampledata.