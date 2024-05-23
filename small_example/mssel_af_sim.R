trajs <- list()
sel.parts <- list()
eff_sizes <- numeric(0)
ss <- numeric(0)
curr.freqs <- numeric(0)


#need a constant that reflects a harmonic series from 1 to 2N-1, but only the
#terms where i/(2N) \in [.01, .99]. for N large this will be very close to the 
#analogous integral, which is ln(.99) - ln(.01) ~= 4.595
ser <- 1:(2*N-1)
harm.const <- sum(1/ser[ser >= .01*2*N & ser <= .99*2*N])


#Simulate and write n.loci allele-frequency trajectories. None of these should have fixed,
#so reject and do it again if fixed. Reject if minor allele frequency is < 0.01.
shift.achieved <- 0
while(shift.achieved == 0){
    pre.sel.freqs <- numeric(0)
    post.sel.freqs    <- numeric(0)
    for(i in 1:n.loci){
        fix <- 1
        while(fix == 1){
            eff_size <- rnorm(1,0,sqrt(herit * sd.trait^2 * harm.const /n.loci))
            s_locus <- eff_size * sel.intens / sd.trait #Charlesworth B3.7.7
            ss[i] <- s_locus
            p0 <- gen.neut.sfs.condit(1, 2*N, .01) #select from neutral sfs conditional on common.
            pre.sel.freqs[i] <- p0 
            sel.part <- sel.traj(p0, s = s_locus, N = N, t = t-t.off)
            post.sel.freqs[i] <- sel.part[nrow(sel.part),2]
            if(t.off > 0){
                post.drift <- neut.traj.time(p0 = sel.part[nrow(sel.part),2], N, t.off)
                sel.part <- rbind(sel.part, cbind(post.drift[,1] + t - t.off, post.drift[,2] ))            
            }
            if(sel.part[nrow(sel.part),2] >= 0.01 & sel.part[nrow(sel.part),2] <= 0.99){fix = 0}
        }
        eff_sizes[i] <- eff_size 
        curr.freqs[i] <- sel.part[nrow(sel.part),2]
        sel.parts[[i]] <- sel.part
    }    
    #Check whether the achieved shift is close to the target.
    trait.sd <- sqrt(sum(2 * curr.freqs * (1 - curr.freqs) *eff_sizes^2))
    sel.shift <- sum(2 * eff_sizes * (post.sel.freqs - pre.sel.freqs)) / trait.sd
    target.shift <- (sel.intens * sd.trait) * (t - t.off) * 2 * N
    if(target.shift*.95 <= sel.shift & sel.shift <= target.shift*1.05){shift.achieved <- 1}
#    print("trait.attempted")
}


for(i in 1:n.loci){
    sel.part <- sel.parts[[i]]
    driftup <- neut.traj(pre.sel.freqs[i], N, loss = TRUE)
    traj.fwd <- rbind(cbind(driftup[,1], rev(driftup[,2])), cbind(sel.part[-1,1] + max(driftup[,1]), sel.part[-1,2] ))
    traj <- traj.fwd
    traj[,2] <- rev(traj.fwd[,2])
    trajs[[i]] <- traj
}

eff_sizes_unscaled <- eff_sizes
eff_sizes <- eff_sizes / trait.sd #The SD of the polygenic score is set
#to be 1 in the present


rm(ser)
rm(harm.const)
rm(eff_size)
rm(s_locus)
rm(p0)
rm(sel.part)
#rm(post.drift)
#rm(post.sel.freqs)
#rm(pre.sel.freqs)


mat.trajs <- matrixify.list.of.trajs(trajs) #matrix of derived allele freq at each locus

phen.traj <- as.numeric(2 * eff_sizes %*% t(mat.trajs))  #phenotypes vector
pt.time <- seq(0, by = 1/(2*N), length.out = max(sapply(trajs, length)/2) ) #desired time interval


n_ders <- numeric(n.loci)
nsites <- numeric(n.loci)
theta_ls <- numeric(n.loci)

ms_trees_list <- list()
ms_haplotypes_list <- list()

print("run mssel")
for(i in 1:n.loci){
    write.traj(traj.fn, trajs[[i]])
    curr.freq.der <- curr.freqs[i]
    n_der <- rbinom(1, n_chroms, curr.freq.der) #number of chroms with derived allele.
    #n_der <- rbinom(1, 2000, curr.freq.der)
    if(n_der == 0){n_der <- n_der + 1}
    if(n_der == n_chroms){n_der <- n_der - 1}
    #if(n_der == 2000){n_der <- n_der - 1}
    n_ders[i] <- n_der
    
    nsite <- 4*N*(1 - dbinom(0,len_hap, re)) 
    nsites[i] <- nsite    
    
    theta <- 4*N* (1 - dbinom(0, len_hap, u))
    theta_ls[i] <- theta    

    ms.string <- paste(ms_dir, "mssel ", as.character(n_chroms), " 1 ", as.character(n_chroms - n_der), " ", as.character(n_der), " ", traj.fn, " ", sel_site, " -r ", as.character(nsite), " ", as.character(len_hap), " -t ", as.character(theta), " -T -L > ", msout.fn, sep = "")
    #ms.string <- paste(ms_dir, "mssel 2000 1 ", as.character(2000 - n_der), " ", as.character(n_der), " ", traj.fn, " ", sel_site, " -r ", as.character(nsite), " ", as.character(len_hap), " -t ", as.character(theta), " -T -L > ", msout.fn, sep = "")
    system(ms.string)

    msoutlines <- readLines(msout.fn)
    dists <- numeric(0)
    for(m in 1:(length(msoutlines) - 4)){
        suppressWarnings(dists[m] <- as.numeric(substring(strsplit(msoutlines[m + 4], "]")[[1]][1],2)))
    }
    dists <- dists[!is.na(dists) & dists != Inf]
    posits.rs <- cumsum(dists)
    ms_to_extract <- which.max(posits.rs[posits.rs < sel_site]) + 1
    if(length(ms_to_extract)  == 0){ms_to_extract <- 1}
    ms_trees_list[[i]] <- read.tree(text = msoutlines[4 + ms_to_extract])
    
    print(paste("simulation locus", as.character(i), "of trait", as.character(phen_nums[k]), "complete"))
}

rm(n_der)
rm(nsite)
rm(theta)

check_samecor_sel_site("ms", ms_trees_list)
split_tree("ms", ms_trees_list)
#store coalesenct times of the whole tree, anc tree and der tree
coal_time_ls(ms_trees_list, anc_trees_ms, der_trees_ms, "ms")
# store the number of lineages of ref allele (col1) and alt allel (col2) into a list of matrix
lins_ls(ms_trees_list, "ms")


#fn_str <- paste("loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_num), sep = "")
fn_str <- paste("loci", as.character(n.loci), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(n_chroms), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(phen_nums[k]), sep = "")

print(fn_str)
save.image(paste(out_dir, "/", fn_str, ".RData", sep = ""))
