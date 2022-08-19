- Unzip the package to any directory.
- The data files are under \data, in particular, syn-lr.mat.
- Add the directories, \inexact_alm_rpca, \lightspeed, \PORPACK, \tensor_toolbox, to your Matlab path.
- Also add \code and its subdirectories to Matlab path.
- To run the experiments on random tensor data, call 
	[Err, Cpu, Iter] = run_batch_tensor_exp( dataname )
	where dataname = 'syn-lr'.  This will produce plots similar to Fig. 4.1 and 4.2 in the paper "Robust Low-rank Tensor Recovery - Models and Algorithms" https://sites.google.com/site/tonyqin/research
	
- The main function for the algorithms is 
	test_trpca( dataname, alg, rRatio, lambdaS, IsTC, verbose, mode, mu1fac, ks ).
	It runs a selected algorithm on a given data set.  Algorithms are indexed.  Refer to get_algs() for the corr. names.  Right now, only alg 1 to 3 are active.  lambdaS = parameter for nuclear norms.  Param for L1 norm = lambdaS*1/sqrt(max_dim)*rRatio.
	IsTC = true for missing data.  mode is for RPCA only, indicating which mode RPCA is applied to.

- syn-lr.mat contains a struct 'data':
data = 

           X: [50x50x20 tensor]
         mag: [1 5]
      obsPct: [1x20 double]
         obs: {1x20 cell}
    noisePct: [0.0100 0.0500 0.1000 0.1500 0.2000 0.2500 0.3000 0.3500 0.4000 0.4500]
       noise: {1x10 cell}
        Nrep: 5

	To get the noisy data with certain % observations and % noise, call gen_syn_data().  You can refer to run_batch_tensor_exp() for calling syntax.

- To run the experiments for non-convex algs, call run_batch_tensor_exp_ncx( dataname ).  The code is similar to the convex case.