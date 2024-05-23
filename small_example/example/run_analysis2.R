library(argparser)
p<- arg_parser("simulation parameters")
p <- add_argument(p, "--n_chromss", help = "number of chromosomes", type = "numeric")
p <- add_argument(p, "--n_loci", help = "number of loci", type = "numeric")
p <- add_argument(p, "--out", help = "out directory", type = "character")
p <- add_argument(p, "--iter_start", help = "start iteration", type = "numeric")
p <- add_argument(p, "--iter_end", help = "end iteration", type = "numeric")
#p <- add_argument(p, "--mssel", help = "run mssel", type = "character")
p <- add_argument(p, "--rent", help = "run rent+", type = "character")
p <- add_argument(p, "--relate", help = "run relate", type = "character")
p <- add_argument(p, "--tsinfer", help = "run tsinfer", type = "character")
p <- add_argument(p, "--argweaver", help = "run argweaver", type = "character")
p <- add_argument(p, "--argneedle", help = "run argneedle", type = "character")
p <- add_argument(p, "--save", help = "save path", type = "character")
analyze_argv <- parse_args(p)

run_mssel <- analyze_argv$mssel
run_rent <- analyze_argv$rent
run_relate <- analyze_argv$relate
run_tsinfer <- analyze_argv$tsinfer
run_aw <- analyze_argv$argweaver
run_an <- analyze_argv$argneedle

#complete phenotype simulations by looping through parameter values and calling 
#pheno_sim_1iter.R for each set of parameters.
N <- 10000
herit <- 1

if(strsplit(analyze_argv$out, "_")[[1]][1] == "nosel"){
    sel.intenses <- .000
    ts <- 0.02
}else if(strsplit(analyze_argv$out, "_")[[1]][1] == "recent"){
    sel.intenses <- .005
    ts <- 0.04
}
n.locis <- analyze_argv$n_loci
n_chromss <- analyze_argv$n_chromss

iter_start <- analyze_argv$iter_start
iter_end <- analyze_argv$iter_end

out <- analyze_argv$out
save <- analyze_argv$save

t.offs <- 0.02
phen_nums <- 1:(iter_end - iter_start + 1) #1 number for each rep we want to do at each combination of parameters

aw_sample <- 50

helper_fn <- "../../helper_functions_coal_sel.R"
#source(helper_fn) #read in helper functions and load ape package

print("load helper")

options(scipen = 999) #disable scientific notation so that parameters passed to ms are given as numbers

time <- c(seq(0, 4, by = 0.001))
pars <- expand.grid(sel.intenses, n.locis, n_chromss, ts, t.offs, phen_nums)

err.array <- array(dim = c(length(time),3,dim(pars)[1]))   
err.std.array <- array(dim = c(length(time),3,dim(pars)[1]))  
qxtest_mat <- matrix(-1, nrow = dim(pars)[1], ncol = 12)


err.array.rent <- array(dim = c(length(time),3,dim(pars)[1]))   
err.std.array.rent <- array(dim = c(length(time),3,dim(pars)[1]))  
qxtest_mat_rent <- matrix(-1, nrow = dim(pars)[1], ncol = 12)


err.array.relate <- array(dim = c(length(time),3,dim(pars)[1]))   
err.std.array.relate <- array(dim = c(length(time),3,dim(pars)[1]))  
qxtest_mat_relate <- matrix(-1, nrow = dim(pars)[1], ncol = 12)


#if(!is.na(analyze_argv$tsinfer)){
	err.array.tsinfer <- array(dim = c(length(time),3,dim(pars)[1]))
	err.std.array.tsinfer <- array(dim = c(length(time),3,dim(pars)[1]))
	qxtest_mat_tsinfer <- matrix(-1, nrow = dim(pars)[1], ncol = 12)
#}

#if(!is.na(analyze_argv$argweaver)){
	err.array.aw <- array(dim = c(length(time), 3, dim(pars)[1]))
	err.std.array.aw <- array(dim = c(length(time),3,dim(pars)[1]))
	qxtest_mat_aw <- matrix(-1, nrow = dim(pars)[1], ncol = 12)
#}


err.array.an <- array(dim = c(length(time), 3, dim(pars)[1]))
err.std.array.an <- array(dim = c(length(time), 3, dim(pars)[1]))
qxtest_mat_aw <- matrix(-1, nrow = dim(pars)[1], ncol = 12)

mat.true.phentrajs <- matrix(nrow = length(time), ncol = 100)

for(iteration in 1:(iter_end - iter_start + 1)){
    fn <- paste(out, "/loci", as.character(n.locis), "_sintens", as.character(round(sel.intenses,3)), "_N", as.character(N), "_nchr", as.character(n_chromss), "_ton", as.character(ts), "_toff", as.character(t.offs), "_herit", as.character(herit), "_", as.character(iteration + iter_start - 1), "_locus1_100.RData", sep = "")    
    print(fn)
    #fn <- paste("recent_20_20_out_full/loci20_sintens0.005_N10000_nchr20_ton0.04_toff0.02_herit1_", as.character(iter), ".RData",  sep = "")
    #fn <- paste("nosel_20_20_out/loci20_sintens0_N10000_nchr2000_ton0.02_toff0.02_herit1_", as.character(iter), ".RData", sep = "")
    load(fn)
    #n_ders <- n_ders[sample_loci]
    #if(!is.na(new_argv$sample_loci_num)){
    #    eff_sizes <- eff_sizes[sample_loci]
    #    eff_sizes <- eff_sizes/0.55
        #trajs <- trajs[sample_loci]
    #    mat.trajs <- matrixify.list.of.trajs(trajs)
    #    phen.traj <- as.numeric(2 * eff_sizes %*% t(mat.trajs)) 
    #    pt.time <- seq(0, by = 1/(2*N), length.out = max(sapply(trajs, length)/2) ) #desired time interval
    #}
    #if(!is.na(run_mssel)){
    source("../analyze_sim_true.R")
    #} 
    if(!is.na(run_rent)){
	source("../analyze_sim_rent.R")
      }
      if(!is.na(run_relate)){
	source("../analyze_sim_relate.R")
      }
      if(!is.na(run_tsinfer)){
	source("../analyze_sim_tsinfer.R")
      }
      if(!is.na(run_aw)){
	print('start to analyze argweaver trees')
	source("../analyze_sim_argweaver.R")
      }
      if(!is.na(run_an)){
	print('start to analyze argneedle trees')
        source("../analyze_sim_argneedle.R")
      }
    
    print(paste("trial", as.character(iteration + iter_start - 1), "complete."))
}
warnings()
#save.image("rent_20_100_analyzed.RData")
save.image(paste(save, "/", strsplit(out, "_")[[1]][1], "_", as.character(n_chromss), "_", as.character(n.locis), "_iter", as.character(iter_start), "_", as.character(iter_end),"_analyzed_trees.RData", sep = ""))
