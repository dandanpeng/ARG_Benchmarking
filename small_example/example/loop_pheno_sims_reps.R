library(ape)
library(reticulate)
library(argparser)
library(data.table)
use_python("/spack/2206/apps/linux-centos7-x86_64_v3/gcc-11.3.0/python-3.11.3-gl2q3yz/bin/python3")
#use_python("/home1/dandanpe/anaconda3/bin/python3.9")
source_python("utils.py")

p<- arg_parser("simulation parameters")
p <- add_argument(p, "--n_chromss", help = "number of chromosomes", type = "numeric")
p <- add_argument(p, "--n_loci", help = "number of loci", type = "numeric")
p <- add_argument(p, "--out", help = "out directory", type = "character")
p <- add_argument(p, "--temp", help = "temp directory", type = "character")
#p <- add_argument(p, "--sample_loci_num", help = "subsample number from mssel output", type = "numeric")
p <- add_argument(p, "--sample_seq_num", help = "subsample sequence numebr from mssel output", type = "numeric")
p <- add_argument(p, "--loci_start", help = "the number of start loci in parrallele running", type = "numeric")
p <- add_argument(p, "--loci_end", help = "the number of end loci in parrallele running", type = "numeric")
p <- add_argument(p, "--rent", help = "run rent+")
p <- add_argument(p, "--relate", help = "run relate")
p <- add_argument(p, "--tsinfer", help = "run tsinfer")
p <- add_argument(p, "--argweaver", help = "run argweaver")
p <- add_argument(p, "--argneedle", help = "run argneedle")
p <- add_argument(p, "--singer", help = "run singer")
#p <- add_argument(p, "--phen_start", help = "phenotype number start", type = "numeric")
#p <- add_argument(p, "--phen_end", help = "phenotype number end", type = "numeric")
p <- add_argument(p, "--iter_start", help = "phenotype number start", type = "numeric")
p <- add_argument(p, "--iter_end", help = "phenotype number end", type = "numeric")
new_argv <- parse_args(p)

#if(!is.na(new_argv$tsinfer)){source_python("ts_py_working.py")}

#complete phenotype simulations by looping through parameter values and calling 
#pheno_sim_1iter.R for each set of parameters.
#N <- 10000
#herit <- 1
out_dir <- new_argv$out
system(paste("mkdir", out_dir))

#sel.intenses <- .005
n.locis <- new_argv$n_loci
n_chromss <- new_argv$n_chromss
sample_seq_num <- new_argv$sample_seq_num
#ts <- 0.04
#t.offs <- 0.02
#phen_nums <- argv$phen_start:argv$phen_end #1 number for each rep we want to do at each combination of parameters


new_temp <- new_argv$temp
system(paste("mkdir", new_temp))
#traj.fn <- paste(temp,"/temp.txt", sep = "")
#msout.fn <- paste(temp1, "/ms_out.txt", sep = "")
#rent_in_fn <- paste(temp, "/rent_in.txt", sep = "")
#relate_haps <- "data.haps"
#relate_sample <- "data.sample"
#relate_map <- "data.map"
#relate_out <- "data"
#argweaver_in_fn <- paste(temp, "/aw.sites", sep = "")
#argweaver_out_fn <- paste(temp, "/aw.out", sep = "")
aw_sample <- 70

helper_fn <- "../../helper_functions_coal_sel.R"
source(helper_fn) #read in helper functions and load ape package

print("load helper")

ms_dir <- "../../msseldir/"
rentplus_fn <- "../../RentPlus.jar"
#len_hap <- 200000 #the length (in base pairs) of the haplotype -- short because 
#we are not worrying about recombination--just need the sel site.
#sel_site <- 100000 #the position of the selected site in the haplotype
#u <- 2e-8 #the neutral mutation rate per base pair/generation
#re <- 2.5e-8 #the recombination rate per base pair/generation
options(scipen = 999) #disable scientific notation so that parameters passed to ms are given as numbers
#sd.trait <- 1

#time <- c(seq(0, 4, by = 0.001))
#pars <- expand.grid(sel.intenses, n.locis, n_chromss, ts, t.offs, phen_nums)

print("start for loop")
#set.seed(100)
#for(k in 1:dim(pars)[1]){
	#sel.intens <- pars[k,1]
	#n.loci <- pars[k,2]
	#n_chroms <- pars[k,3]
	#t <- pars[k,4]
	#t.off <- pars[k,5]
	#phen_num <- pars[k,6]
	#source("../pheno_sim_1iter.R")
	#$source("../mssel_af_sim.R")
	#print(paste("trial", as.character(phen_num), "complete."))
#}


#save.image(paste("sim_trees", ".RData", sep = ""))




#err.array <- array(dim = c(length(time),3,dim(pars)[1]))   
#err.std.array <- array(dim = c(length(time),3,dim(pars)[1]))  
#qxtest_mat <- matrix(-1, nrow = dim(pars)[1], ncol = 12)

#err.array.rent <- array(dim = c(length(time),3,dim(pars)[1]))   
#err.std.array.rent <- array(dim = c(length(time),3,dim(pars)[1]))  
#qxtest_mat_rent <- matrix(-1, nrow = dim(pars)[1], ncol = 12)

#err.array.relate <- array(dim = c(length(time),3,dim(pars)[1]))   
#err.std.array.relate <- array(dim = c(length(time),3,dim(pars)[1]))  
#qxtest_mat_relate <- matrix(-1, nrow = dim(pars)[1], ncol = 12)

#err.array.tsinfer <- array(dim = c(length(time),3,dim(pars)[1]))   
#err.std.array.tsinfer <- array(dim = c(length(time),3,dim(pars)[1]))  
#qxtest_mat_tsinfer <- matrix(-1, nrow = dim(pars)[1], ncol = 12)


#err.array.aw <- array(dim = c(length(time), 3, dim(pars)[1]))
#err.std.array.aw <- array(dim = c(length(time),3,dim(pars)[1]))
#qxtest_mat_aw <- matrix(-1, nrow = dim(pars)[1], ncol = 12)

#mat.true.phentrajs <- matrix(nrow = length(time), ncol = dim(pars)[1])

for(iter in new_argv$iter_start:new_argv$iter_end){
    print(strsplit(new_argv$out, '_')[[1]][1])
    if(strsplit(new_argv$out, '_')[[1]][1] == "nosel"){
        ms_fn <- paste("nosel_true_sim/2000sample/loci100_sintens0_N10000_nchr2000_ton0.02_toff0.02_herit1_", as.character(iter), ".RData", sep = "")
    }else if(strsplit(new_argv$out, '_')[[1]][1] == "recent"){
        ms_fn <- paste("recent_true_sim/loci100_sintens0.005_N10000_nchr2000_ton0.04_toff0.02_herit1_", as.character(iter), ".RData", sep = "")
    }
    load(ms_fn)
    n.loci <- new_argv$n_loci
    n_chromss <- new_argv$n_chromss
    traj.fn <- paste(new_temp,"/temp.txt", sep = "")
    msout.fn <- paste(new_temp, "/ms_out.txt", sep = "")
    rent_in_fn <- paste(new_temp, "/rent_in.txt", sep = "")
    rent_sample_in <- paste(new_temp, "/rent_sample_in.txt", sep = "")
    argweaver_in_fn <- paste(new_temp, "/aw.sites", sep = "")
    argweaver_out_fn <- paste(new_temp, "/aw.out", sep = "")
    #out_dir <- "recent_20_20_out"
    aw_sample <- 70
    running_time <- c()
    source(helper_fn)
    #if(!is.na(new_argv$argneedle)){
	#asmc_freq(n_ders)
        #create_decoding_quantities("asmc_decoding/NE10K.demo",
        #                           "asmc_decoding/decoding.disc",
        #                           "asmc_decoding/decoding.frq")
    #}
    source("../sample_and_run.R")
    #source("../sample_run_wo_rescale.R")
}

save(running_time, file = "relate_200_run_time.RData")
print(paste("running time:", as.character(sum(running_time))))
#for(iter in 1:dim(pars)[1]){
	#sel.intens <- pars[iter,1]
	#n.loci <- pars[iter,2]
	#n_chroms <- pars[iter,3]
	#t <- pars[iter,4]
	#t.off <- pars[iter,5]
	#phen_num <- pars[iter,6]
	#fn <- paste(out_dir, "/loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), ".RData", sep = "")
	#load(fn)
	#n.loci <- 20
	#n_chroms <- 20
        # source("../sample_and_run.R")
	#source("../analyze_sim_true.R")
	#source("../analyze_sim_rent.R")
	#source("../analyze_sim_relate.R")
	#source("../analyze_sim_tsinfer.R")
	#print(paste("trial", as.character(iter), "complete."))
#}


#save.image(paste("analyzed_trees", ".RData", sep = ""))

#source("make_figures.R")
