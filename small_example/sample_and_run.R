# if(!is.na(run_argv$sample_loci_num)){
#     print(paste("We are going to sample ", as.character(run_argv$sample_loci_num), "loci."))
#     sample_loci <- sort(sample(1:100, run_argv$sample_loci_num))
#     if(sel.intenses != 0){
#     	sample.pre.sel.freqs <- pre.sel.freqs[sample_loci]
#     	sample.post.sel.freqs <- post.sel.freqs[sample_loci]
#     
#     	sample.eff_sizes <- eff_sizes[sample_loci]
#     	sample.curr.freqs <- curr.freqs[sample_loci]
#     
#     	sample.trait.sd <- sqrt(sum(2 * sample.curr.freqs * (1 - sample.curr.freqs) * sample.eff_sizes^2))
#     	sample.sel.shift <- sum(2 * sample.eff_sizes * (sample.post.sel.freqs - sample.pre.sel.freqs)) / sample.trait.sd
#     
#     	# check if the sampled loci are representative enough
#     	while(target.shift * .95 >= sel.shift & sel.shift >= target.shift * 1.05){
#         	sample_loci <- sort(sample(1:100, sample_loci_num))
#         	sample.pre.sel.freqs <- pre.sel.freqs[sample_loci]
#         	sample.post.sel.freqs <- post.sel.freqs[sample_loci]
#         	sample.eff_sizes <- eff_sizes[sample_loci]
#         	sample.curr.freqs <- curr.freqs[sample_loci]
#         	sample.trait.sd <- sqrt(sum(2 * sample.curr.freqs * (1 - sample.curr.freqs) * sample.eff_sizes^2))
#         	sample.sel.shift <- sum(2 * sample.eff_sizes * (sample.post.sel.freqs - sample.pre.sel.freqs)) / sample.trait.sd    
#     	}
#     }	
#     print("Found representative loci samples.")
#     
#     trajs <- trajs[sample_loci]
#     n_ders <- n_ders[sample_loci]
#     nsites <- nsites[sample_loci]
#     theta_ls <- theta_ls[sample_loci]
# }


sample_indices_ls <- list()

ms_trees_list <- list()
rent_trees_list <- list()
relate_trees_list <- list()
tsinfer_trees_list <- list()
argweaver_samples_list <- list()
argweaver_trees_list <- list()
argneedle_trees_list <- list()
singer_samples_list <- list()
singer_trees_list <- list()
ms_haplotypes_list <- list()

for(i in run_argv$loci_start:run_argv$loci_end){
    write.traj(traj.fn, trajs[[i]])
    n_der <- n_ders[i]
    nsite <- nsites[i]
    theta <- theta_ls[i]
    
    #counter = 0
    success = 1
    
    while (success == 1){
        #ms.string <- paste(ms_dir, "mssel ", as.character(2000), " 1 ", as.character(2000 - n_der), " ", as.character(n_der), " ", traj.fn, " ", sel_site, " -r ", as.character(nsite), " ", as.character(len_hap), " -t ", as.character(theta), " -T -L > ", msout.fn, sep = "")
        
        ms.string <- paste(ms_dir, "mssel ", as.character(n_chroms), " 1 ", as.character(n_chroms - n_der), " ", as.character(n_der), " ", traj.fn, " ", sel_site, " -r ", as.character(nsite), " ", as.character(len_hap), " -t ", as.character(theta), " -T -L > ", msout.fn, sep = "")
        system(ms.string)
    
        #cut down mssel output to just the part rent+ needs
        process.string <- paste("sed '/positions: /,$!d' ", msout.fn, " | sed 's/.*positions: //' > ", rent_in_fn, sep = "")
        system(process.string)
        
        
        #################################### RENT+ #################################
        # sample from mssel output (actually transferred mssel ouput which can be fed to RENT+) 
        if(!is.na(run_argv$sample_seq_num)){
            sample_rent_input(rent_in_fn, rent_sample_in, run_argv$sample_seq_num, i - run_argv$loci_start + 1) #sample chromosomes
            
            if(!is.na(run_argv$rent)){
                rent.string <- paste("java -jar ", rentplus_fn, " -t -l ", as.character(len_hap), " ", rent_sample_in, sep = "")
                start.time <- Sys.time()
                success = system(rent.string)
                end.time <- Sys.time()
                running_time <- c(running_time, end.time - start.time)
                if(success == 0){rent_in_fn = rent_sample_in} #when we need to sample and run RENT+
            }else{
                rent_in_fn = rent_sample_in
                success = 0 # when we need to sample but don't want to run RENT+
            }
        }else{ 
            if(!is.na(run_argv$rent)){# when we don't need to sample but need to run RENT+
		start.time <- Sys.time()
                rent.string <- paste("java -jar ", rentplus_fn, " -t -l ", as.character(len_hap), " ", rent_in_fn, sep = "")
                success = system(rent.string)
		end.time <- Sys.time()
		print(paste("running time: ", as.character(end.time - start.time)))
            }else{ #when we don't need to sample and don't need to run RENT+
                success = 0
            }
        }
    }  
      
    
  #pull trees from ms
	msoutlines <- readLines(msout.fn)
	dists <- numeric(0)
	for(m in 1:(length(msoutlines) - 4)){
		suppressWarnings(dists[m] <- as.numeric(substring(strsplit(msoutlines[m + 4], "]")[[1]][1],2)))
	}
	dists <- dists[!is.na(dists) & dists != Inf]
	posits.rs <- cumsum(dists)
	ms_to_extract <- which.max(posits.rs[posits.rs < sel_site]) + 1
	if(length(ms_to_extract)  == 0){ms_to_extract <- 1}
	if(!is.na(run_argv$sample_seq_num)){
	    pre_tree <- read.tree(text = msoutlines[4 + ms_to_extract])
	    ms_trees_list[[i - run_argv$loci_start + 1]] <- drop.tip(pre_tree, setdiff(as.character(c(1:2000)), as.character(sample_indices_ls[[i - run_argv$loci_start + 1]])))
	}else{
	    ms_trees_list[[i - run_argv$loci_start + 1]] <- read.tree(text = msoutlines[4 + ms_to_extract])
	}
    
    if(!is.na(run_argv$rent)){
        print("start to run rent+")
        #Read in trees at selected site in rent+
        rent_trees_fn <- paste(rent_in_fn, ".trees", sep = "")
        rent_trees <- readLines(rent_trees_fn)
    
        #which tree is for the selected site?
        rent_dists <- numeric(0)
        for(r in 1:length(rent_trees)){
            rent_dists[r] <- strsplit(rent_trees[r], "\t")[[1]][1]
        }
        rent_extract_index <- which.min(abs(as.integer(rent_dists) - sel_site))
        rent_to_extract <- strsplit(rent_trees[rent_extract_index], "\t")[[1]][1]    

        rent_sel_ind <- grep(paste(rent_to_extract, "\t", sep = ""), rent_trees)
    
        if(length(rent_sel_ind) > 1){
            rent_sel_ind <- rent_sel_ind[1]
        }
    
        rent_sel_tree_newick <- sub(paste(rent_to_extract, "\t", sep = ""), "", rent_trees[rent_sel_ind])
        rent_trees_list[[i - run_argv$loci_start + 1]] <- read.tree(text = paste(rent_sel_tree_newick, ";", sep = ""))
    }
    
    snp_pos <- recover_pos(rent_in_fn, len_hap)
    sel_site <- snp_pos[which.min(abs(snp_pos - sel_site))]
    print(sel_site)
    derived_info <- hap_to_mat(rent_in_fn)

    ############################### RELATE #############################
    #transfer rent+ input to RELATE recognized format
    if(!is.na(run_argv$relate)){
        relate_input(derived_info, snp_pos, n_chromss)
        
        print("start to run RELATE")
       
        setwd(new_temp)
	options(digits.secs = 6)
	start.time <- Sys.time()
        relate.string <- paste("../relate/Relate --mode All -m", as.character(u),
                               "-N", as.character(N*2),
                               "--haps data.haps",
                               "--sample data.sample",
                               "--map data.map",
                               "--output data")
        system(relate.string)
        end.time <- Sys.time()
	#print(paste("running time", as.character(end.time - start.time)))
        #Read in trees at selected site in RELATE
        relate_sel_tree <- paste("../relate/RelateExtract --mode AncToNewick",
                              "--anc data.anc",
                              "--mut data.mut",
                              "--bp_of_interest", as.character(sel_site),
                              "--output data")
        print("execute relateextract")
        system(relate_sel_tree)
        relate_trees_list[[i - run_argv$loci_start + 1]] <- read.tree(text = readLines("data.newick"))
        setwd("..")
        
        for(k in 1:length(relate_trees_list[[i - run_argv$loci_start + 1]]$tip.label)){
          relate_trees_list[[i - run_argv$loci_start + 1]]$tip.label[k] <- as.character(as.numeric(relate_trees_list[[i - run_argv$loci_start + 1]]$tip.label[k]) + 1)
        }
    }
    
    ################################ tsinfer #############################
    if(!is.na(run_argv$tsinfer)){
        #run tsinfer
        print("start to run tsinfer")
        #source("TSinfer.R")
        data =  cbind(snp_pos, derived_info)
	      #write.csv(data, paste(run_argv$temp, "/data.csv", sep = ''))
	      #sel_row = which(data[, 1] == sel_site)

        #Read in tree at selected site in tsinfer
	options(digits.secs = 6)
	start.time <- Sys.time()
	      tsinfer_sel_tree_newick <- tsfunc(data, sel_site - 1, N, u) # sel_site - 1 because the position is 0-based in tsinfer 
	end.time <- Sys.time()
	print(paste("running time", as.character(end.time - start.time)))
        tsinfer_trees_list[[i- run_argv$loci_start + 1]] <- read.tree(text = tsinfer_sel_tree_newick)
      	tsinfer_trees_list[[i - run_argv$loci_start + 1]]$tip.label <- as.character(as.numeric(gsub("n", "", tsinfer_trees_list[[i - run_argv$loci_start + 1]]$tip.label)) + 1)
    }
    
    ###############################  ARGWeaver #############################
  
    if(!is.na(run_argv$argweaver)){
        print("start to run ARGWeaver")
        
        aw_input(rent_in_fn, snp_pos, len_hap, new_temp) 
        
        #run ARGweaver
	start.time <- Sys.time()
        argweaver.string <- paste("/home1/dandanpe/sc1/argweaver_new/bin/arg-sample", "-s", argweaver_in_fn, "-N", N, "-r", re, "-m", u, "-n", as.character(aw_sample * 10), "--smc-prime --overwrite -o", argweaver_out_fn)
        system(argweaver.string)
        end.time <- Sys.time()
	print(end.time - start.time)
        #Read in tree at selected site in ARGweaver
        for(sample in 0:(aw_sample-20)){
	    real_sample <- sample + 20
            aw.smc <- read.table(paste(argweaver_out_fn, ".", as.character(10*real_sample), ".smc.gz", sep = ""), header = FALSE, fill = TRUE)
            aw.smc.trees <- aw.smc[3:nrow(aw.smc), ]
            aw.smc.trees$V2 <- as.numeric(aw.smc.trees$V2)
            aw.smc.trees$V3 <- as.numeric(aw.smc.trees$V3)
            sel_row <- aw.smc.trees[aw.smc.trees$V1 == "TREE" & aw.smc.trees$V2 <= sel_site & aw.smc.trees$V3 >= sel_site, ]
            sel_tree <- read.tree(text = sel_row$V4)
            real_label <- as.vector(unlist(aw.smc[1, 2:length(aw.smc)]))
            sel_tree$tip.label <- unlist(lapply(sel_tree$tip.label, FUN = change_aw_label))
            argweaver_samples_list[[sample + 1]] <- sel_tree
        }
        argweaver_trees_list[[i - run_argv$loci_start + 1]] <- argweaver_samples_list
    }
    
    ###############################  ARG-Needle #############################
    
    if(!is.na(run_argv$argneedle)){
      print("start to run ARG-Needle")
      arg_needle_input(derived_info, snp_pos, n_chromss)
      asmc_freq(derived_info, snp_pos)
      create_decoding_quantities("asmc_decoding/NE10K.demo", 
                                 "asmc_decoding/decoding.disc",
                                 "asmc_decoding/decoding.frq")
      setwd(new_temp)
      start.time <- Sys.time()
      system("~/.local/bin/arg_needle --asmc_clust 1 --normalize_demography ../asmc_decoding/NE10K.demo --hap_gz an_data.hap --map an_data.map --asmc_decoding_file custom.decodingQuantities.gz --out an_data")
      end.time<- Sys.time()
      print(end.time - start.time)
      setwd("..")
      an_newick(paste(new_temp, "/an_data.argn", sep = ""), paste(new_temp, "/argn_newick.txt", sep = ""))
      argn_tree = read.tree(file = paste(new_temp, "/argn_newick.txt", sep = ''))
      argn_tree$tip.label <- as.character(as.numeric(argn_tree$tip.label) + 1)
      argneedle_trees_list[[i - run_argv$loci_start + 1]] <- argn_tree
      
    }
    
    ###############################  Singer #############################
    if(!is.na(run_argv$singer)){
      print("start to run Singer")
      
      setwd(new_temp)
      vcf_input(derived_info, snp_pos)

      start.time <- Sys.time()
      system(paste("~/proj/singer-0.1.7-beta-linux-x86_64/singer_master -Ne 10000 -m 2e-8 -ratio 1.25 -start 0 -end", as.character(len_hap-1), "-vcf singer -output singer -n 700 -thin 10 -polar 0.99"))
      end.time <- Sys.time()
      print(end.time - start.time)
      system("~/proj/singer-0.1.7-beta-linux-x86_64/convert_to_tskit -input singer -output singer -start 200 -end 699 -step 10")
      for(index in 1:50){
        tree_index <- (19 + index) * 10
        treeSeqFile <- paste("singer_", as.character(tree_index), ".trees", sep = "")
        singer_samples_list[[index]] <- read.tree(text = get_singer_tree(treeSeqFile, sel_site))
        singer_samples_list[[index]]$tip.label <-  as.character(as.numeric(gsub("n", "", singer_samples_list[[index]]$tip.label)) + 1)
      }
      singer_trees_list[[i - run_argv$loci_start + 1]] <- singer_samples_list

      
      setwd("..")
    }
    
    #save the haplotypes produced by ms (and fed to rent+) in a list entry.
    #the result will be a character vector. The first entry in the vector is relative (0-1) positions, with the selected
    #site at 0.05. The remaining entries are string haplotypes of 0s and 1s.ts
    ms_haplotypes_list[[i - run_argv$loci_start + 1]] <- readLines(rent_in_fn)
    
    print(paste("simulation locus", as.character(i), "of trait", as.character(iter), "complete"))
}

save.image("check_trees.RData")

##rescale branch length of RENT+ trees
if(!is.na(run_argv$rent)){
  for(i in 1:(run_argv$loci_end - run_argv$loci_start + 1)){
    rent_trees_list[[i]]$edge.length <- rent_trees_list[[i]]$edge.length*2
  }
}

##rescale branch length of relate trees
if(!is.na(run_argv$relate)){
    for(i in 1:(run_argv$loci_end - run_argv$loci_start + 1)){
      relate_trees_list[[i]]$edge.length <- relate_trees_list[[i]]$edge.length/(28 * 2 * N)
    }
}


##rescale branch length of tsinfer trees
##The original branch lengths are given in units of generations.
if(!is.na(run_argv$tsinfer)){
  for(i in 1:(run_argv$loci_end - run_argv$loci_start + 1)){
      tsinfer_trees_list[[i]]$edge.length <- tsinfer_trees_list[[i]]$edge.length*2
  }
}

##rescale branch length of argweaver trees
##The original branch lengths are given in units of generations.
if(!is.na(run_argv$argweaver)){
    for(i in 1:(run_argv$loci_end - run_argv$loci_start + 1)){
        for(j in 1:length(argweaver_trees_list[[i]])){
            argweaver_trees_list[[i]][[j]]$edge.length <- argweaver_trees_list[[i]][[j]]$edge.length/(2*N)
        }
    }
}

##rescale branch length of argneedle trees
if(!is.na(run_argv$argneedle)){
  for(i in 1:(run_argv$loci_end - run_argv$loci_start + 1)){
    argneedle_trees_list[[i]]$edge.length <- argneedle_trees_list[[i]]$edge.length/(2*N)
  }
}

##rescale branch length of singer trees
if(!is.na(run_argv$singer)){
  for(i in 1:(run_argv$loci_end - run_argv$loci_start + 1)){
    for(j in 1:length(singer_trees_list[[i]])){
      singer_trees_list[[i]][[j]]$edge.length <- singer_trees_list[[i]][[j]]$edge.length/(2*N)
    }
  }
}



#Run through list and check for sites where 2+ variants had same coordinates as
#selected site. At those places, check which tree(s) are monophyletic for derived tips.
#if more than one, select randomly from among them.
#split into ancestral and derived trees and retrieve coalescence times.
check_samecor_sel_site("ms", ms_trees_list)
split_tree("ms", ms_trees_list)
#store coalesenct times of the whole tree, anc tree and der tree
coal_time_ls(ms_trees_list, anc_trees_ms, der_trees_ms, "ms")
# store the number of lineages of ref allele (col1) and alt allel (col2) into a list of matrix
lins_ls(ms_trees_list, "ms")


if(!is.na(run_argv$rent)){
    check_samecor_sel_site("rent", rent_trees_list)
    split_tree("rent", rent_trees_list)
    coal_time_ls(rent_trees_list, anc_trees_rent, der_trees_rent, "rent")
    lins_ls(rent_trees_list, "rent")
}


if(!is.na(run_argv$relate)){
    check_samecor_sel_site("relate", relate_trees_list)
    split_tree("relate", relate_trees_list)
    coal_time_ls(relate_trees_list, anc_trees_relate, der_trees_relate, "relate")
    lins_ls(relate_trees_list, "relate")
}

if(!is.na(run_argv$tsinfer)){
    check_samecor_sel_site("tsinfer", tsinfer_trees_list)
    split_tree("tsinfer", tsinfer_trees_list)
    coal_time_ls(tsinfer_trees_list, anc_trees_tsinfer, der_trees_tsinfer, "tsinfer")
    lins_ls(tsinfer_trees_list, "tsinfer")
}

if(!is.na(run_argv$argweaver)){
    mcmc_check_samecor_sel_site(argweaver_trees_list)
    mcmc_split_tree("argweaver", argweaver_trees_list)
    mcmc_coal_time_ls(argweaver_trees_list, anc_trees_argweaver, der_trees_argweaver, "argweaver")
    aw_lins_ls(argweaver_trees_list, "argweaver")
}

if(!is.na(run_argv$argneedle)){
  check_samecor_sel_site("argneedle", argneedle_trees_list)
  split_tree("argneedle", argneedle_trees_list)
  coal_time_ls(argneedle_trees_list, anc_trees_argneedle, der_trees_argneedle, "argneedle")
  lins_ls(argneedle_trees_list, "argneedle")
}

if(!is.na(run_argv$singer)){
  mcmc_check_samecor_sel_site(singer_trees_list)
  mcmc_split_tree("singer", singer_trees_list)
  mcmc_coal_time_ls(singer_trees_list, anc_trees_singer, der_trees_singer, "singer")
  mcmc_lins_ls(singer_trees_list, "singer")
}


fn_str <- paste("loci", as.character(length(trajs)), "_sintens", as.character(round(sel.intens,3)), "_N", as.character(N), "_nchr", as.character(run_argv$n_chromss), "_ton", as.character(t), "_toff", as.character(t.off), "_herit", as.character(herit), "_", as.character(iter), "_locus", as.character(run_argv$loci_start), "_", as.character(run_argv$loci_end), sep = "")


save.image(paste(run_argv$out, '/', fn_str, ".RData", sep = ""))


	


    

