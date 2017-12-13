```
This application estimates the selfing rates for subpopulations and
 classfies each individual into subpopulations given the sequence data.
 
Synopsis:
	inbreed -d data_file -o output_file [-i initial_file] 
	[-K population number] [-L loci number] [-N total individual number]
	[-p ploid] [-u iteration number] [-b burn-in number] [-m missingdata] 
	[-t thinning] [-c chain number] [-s seed1 seed2 seed3] [-sl significance level] 
	[-b label] [-a popdata] [-g GR_flag] [-r ckrep] [-f prior_flag] [-v pop_ind_flag]

Parameters:
-d data_file - name of data file
-o output_file - name of output file
-i initial_file - name of initial file
-K popnum - subpopulation number
-L locinum - totoal loci number
-N totalsize - total individual number
-p ploid -  the number of haplotype in a genome
-u update - MCMC iteration number
-b burnin - MCMC burn-in number
-t thinning - MCMC thinning interval length
-c chainnum - the number of MCMC chains
-s seed1 seed2 seed3 three integers to override the default seed selection 
-m an integer represents missing data 
-sl significant level for the confidence intervals in result_analysis.c
-b label boolean indicates whether data_file contains labels for individuals, 1=yes, 0=no
-a popdata boolean indicates whether data_file contains a column about the original population information, 1=yes, 0=no
-g GR_flag boolean indicates whether Gelman_Rudin statistic is used to check convergence,1=yes, 0=no
-r ckrep integer indicates how many stored iterations after burn-in are used in convergence checking
-f prior_flag boolean indicates which prior for selfing rates, uniform (0) or normal (1) or DPM (2) Dirichlet Process prior
-v pop_ind_flag boolean indicates whether selfing rates are for populations (0) or individuals (1)
 
In file inbreed.c:
	main() :
	calls the following functions:
		param_decomp() :decomposition of the commandline arguments
		read_seqs() in data_interface.c:opens the datafile and reads in the sequence data
		read_init() in initial.c:read in the initial values of inbreeding coefficients (F) for MCMC updating
		printinfo() :print the basic information of the running condition
		mcmc_updating_POP() in mcmc.c:implements MCMC updating for admixture models and estimates the selfing rates of subpopulations
		mcmc_updating_INDV() in mcmc.c:implements MCMC updating for admixture models and estimates the selfing rates of individuals
		chain_converg() in check_converg.c:check the convergence of selfing rates using the Gelman_Rubin statistics and print the mixing condition to the output file
		chain_stat() in result_analysis.c:summarize all the results from MCMC chains and print the final results to output file
	printinfo() :print the basic information of the running condition	
	param_decomp() :decomposition of the commandline arguments


In file data_interface.c:
	read_seqs() :opens the datafile and reads in the sequence data
	calls the following functions:
		isnew();
		word_cnt();
		word_split();
		get_missing();
	isnew() :checks whether array has the element of the same
		 value as "num"
	word_cnt() :counts the number of words in a string
	word_split() :splits the string into "num" words and makes 
		      each word an integer number
	get_missing() :generate the matrix and vector that stores the missing data

In file initial.c:
	read_init() :read in the initial values of inbreeding coefficients 
		     (F) for MCMC updating
	calls the following functions:
		word_cnt();
		word_split();
		int_split();
		int_to_string();
	word_cnt() :counts the number of words in a string
	word_split() :splits the string into "num" words and makes 
		      each word a floating number
	int_split()
	int_to_string() :takes in an integer number and change it into a string for print

In file mcmc.c:
	mcmc_updating_POP() :implements MCMC updating for admixture models
	calls the following functions:
		allocate_chn();
		allocate_space();
		update_P();
		update_SG_POP();
		update_ZQ();
		update_alpha();
		move_nodes();
		print_info();
		store_chn();
	mcmc_updating_INDV() :implements MCMC updating for admixture models
	calls the following functions:
		allocate_chn();
		allocate_space();
		init_SG_dp();
		update_P();
		update_S_dp();
		update_S();
		update_G();
		update_ZQ();
		update_alpha();
		move_nodes();
		print_info();
		store_chn();
	allocate_space() :Allocate the space for UPMCMC structure
	allocate_chn() : Allocate the space for CHAINS structure
	init_SG_dp() :Initially draw the selfing rate of each individual from the Dirichlet Process Prior
	calls the following functions:
		insert();
		creat();
		delete();
		disc_unif();
		rgeom();
	update_P() :updating P (allele frequency per locus per subpop) with Dirichlet distribution
	calls the following functions:
		max();
		rdirich();
	update_S() :Update the selfing rates of individuals using either uniform or normal prior
	calls the following functions:
		sample_mu_sigma();
		dgeom();
		loglkd_s2();
	update_SG_POP() :Update the selfing rates of subpopulations and store them in "ptr"
	calls the following functions:
		proposal();
		log_ld_indv();
	update_S_dp() :Update the selfing rates with the Dirichlet Process Prior
	calls the following functions:
		insert();
		creat();
		delete();
		disc_unif();
		rgeom();
	update_G() :Update G with the individual selfing rates according to the Geometric distribution and the Metropolis-Hastings ratio
	calls the following functions:
		log_ld_indv();
	update_ZQ() :Update Z and Q according to the Metropolis-Hastings ratio
	calls the following functions:
		disc_unif();
		rdirich();
	update_alpha() :Update alpha according to the Metropolis-Hastings ratio
	move_nodes() :copy the updated information from "ptr" to "tmptr" and store there for the next updating step
	print_info() :print the information of each updating step
	store_chn() :store the updated information in CHAINS at each iteration
	max() :return the largest value in the vector "vec"
	find() :find the node with element "selfing_rate" less than or equal "self" and return the found node address
	creat() :creat a new node and insert it into the list
	calls the following functions:
		find();
	insert() :insert an individual into one of existing nodes
	delete() :delete the individual and any node that its element "num" is zero
	dgeom() :calculate the density value of geometric distribution
	sample_mu_sigma() :draw mu and var from a hierarchical normal prior, mu ~ Normal and var ~ Inverse Gamma
	loglkd_s2() :calculate the log normal prior of selfing rate given the generation
	genofreq() :calculates genotype frequency given selfing 
		    generations and current allele frequencies
	proposal() :calculates the probability of inbreeding coefficients
		    for each population given generations
	chcksame() :checks whether all the elements in the integer vector
		    "pop" of length "num" are the same or not
	chckdbsm() :checks whether all the elements in the double vector
		    "pop" of length "num" are the same or not
	log_ld_indv() :Calculate the log likelihood of a given generation for individual "index"

In file check_converg.c:
	chain_converg() :check the convergence of inbreeding coefficients
			 using the Gelman_Rubin statistics and print the
			 mixing condition to the output file
	calls the following functions:
		GelmanRubin();
		double_comp();
	GelmanRubin() :calculate the the Gelman_Rubin statistics
	double_comp() :compare two double numbers in order to sort the
		       array of doubles.It is called in function qsort()

In file result_analysis.c:
	chain_stat() :summarize all the results from MCMC chains
		      and print the final results to output file
	calls the following functions:
		summstat();
	summstat() :calculate the summary statistics of the number
		    in the vector "vec" of "length"
	calls the following functions:
		quantile() in quantile.c;

In file quantile.c:
	quantile() :seeks the number with certain quantile in an array
	calls the following functions:
		indexx();
	indexx() :sorts the array "arr" and stores the order in "indx"
	calls the following functions:
		ivector();
		free_ivector();
		nrerror();
	ivector() :allocate an int vector with subscript range v[nl..nh]
	free_ivector() :free an int vector allocated with ivector()
	nrerror() :Numerical Recipes standard error handler

In file random.c:
	nrerror() : Numerical Recipes standard error handler
	wichmann() :generate the uniform random variables on (0,1)
	ran1() : generate the uniform random variables on (0,1)
	setseeds() :set the seeds of the random generator
	printseeds() :print the seeds into a file
	gammln() :returns ln[Gamma[xx]] for xx > 0
	factln() :returns ln(n!)
	bico() :returns binomial coefficients (nCk)
	moment() :Calculate the mean and variance of "data
	rexp() :generates from an exponential distribution
	poidev() :Returns a poisson RV with mean xm
	rgamma1() :generates from a gamma distribution with alpha < 1
	rgamma2() :generates from a gamma distribution with alpha > 1
	rgamma() :generates from a general gamma(alpha,beta) distribution
	rbeta() :Generates from a beta (alpha,beta) distribution
	rdirich() :Generates from a Dirichlet distribution
	rstd_normal() :generates from a standard normal(0,1) distribution
	rnormal() :generates from a general normal(mu,sigma) distribution
	rgeom() :generates from a geometric distribution
	rbinom() :Generates from a binomial (N,p) distribution
	bnldev() :
	disc_unif() :generate the integer that represent the interval where a uniform random value falls
```
