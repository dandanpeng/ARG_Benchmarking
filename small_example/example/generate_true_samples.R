library(argparser)

p<- arg_parser("simulation parameters")
p <- add_argument(p, "--n_chromss", help = "number of chromosomes", type = "numeric")
p <- add_argument(p, "--n_loci", help = "number of loci", type = "numeric")
p <- add_argument(p, "--out", help = "out directory", type = "character")
p <- add_argument(p, "--temp", help = "temp directory", type = "character")
#p <- add_argument(p, "--phen_start", help = "phenotype number start", type = "numeric")
#p <- add_argument(p, "--phen_end", help = "phenotype number end", type = "numeric")
p <- add_argument(p, "--iter_start", help = "phenotype number start", type = "numeric")
p <- add_argument(p, "--iter_end", help = "phenotype number end", type = "numeric")
argv <- parse_args(p)

new_argv <- list()
new_argv$loci_end <- 100
new_argv$loci_start <- 1
sample_seq_num <- NA

N <- 10000
herit <- 1
out_dir <- argv$out
system(paste("mkdir", out_dir))

sel.intenses <- 0.005
n.locis <- argv$n_loci
n_chromss <- argv$n_chromss
ts <- 0.04
t.offs <- 0.02
phen_nums <- argv$iter_start:argv$iter_end #1 number for each rep we want to do at each combination of parameters
print(phen_nums)
#sample_num <- 20

temp <- argv$temp
system(paste("mkdir", temp))
traj.fn <- paste(temp,"/temp.txt", sep = "")
msout.fn <- paste(temp, "/ms_out.txt", sep = "")

helper_fn <- "../../helper_functions_coal_sel.R"
source(helper_fn) #read in helper functions and load ape package

print("load helper")

ms_dir <- "../../msseldir/"

len_hap <- 200000 #the length (in base pairs) of the haplotype -- short because 
#we are not worrying about recombination--just need the sel site.
sel_site <- 100000 #the position of the selected site in the haplotype
u <- 2e-8 #the neutral mutation rate per base pair/generation
re <- 2.5e-8 #the recombination rate per base pair/generation
options(scipen = 999) #disable scientific notation so that parameters passed to ms are given as numbers
sd.trait <- 1

time <- c(seq(0, 4, by = 0.001))
pars <- expand.grid(sel.intenses, n.locis, n_chromss, ts, t.offs, phen_nums)
#pars <- expand.grid(sel.intenses, 20, 20, ts, t.offs, 1:30)

for(k in 1:length(phen_nums)){
    print(phen_nums[k])
    sel.intens <- pars[k,1]
    n.loci <- pars[k,2]
    n_chroms <- pars[k,3]
    t <- pars[k,4]
    t.off <- pars[k,5]
    phen_num <- pars[k,6]
    source("../mssel_af_sim.R")
    print(paste("trial", as.character(phen_num), "complete."))
}

