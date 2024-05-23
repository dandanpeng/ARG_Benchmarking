#neutral
trajs_neut_aw <- array(dim = c(length(time), length(argweaver_trees_list), length(argweaver_trees_list[[1]]))) #matrix of proportion of derived lineages (1st estimator)
vars_neut_bin_aw <- array(dim = c(length(time), length(argweaver_trees_list), length(argweaver_trees_list[[1]])))

avg_trajs_neut_aw <- matrix(nrow = length(time), ncol = length(argweaver_trees_list))
avg_vars_neut_bin_aw <- matrix(nrow = length(time), ncol = length(argweaver_trees_list))

for(i in 1:length(argweaver_trees_list)){	
    for(j in 21:length(argweaver_trees_list[[i]])){
    	trajs_neut_aw[,i,j] <- est_af_traj_neut(lins.list.argweaver[[i]][[j]])
    	vars_neut_bin_aw[,i,j] <- est_af_var_neut_bin(lins.list.argweaver[[i]][[j]])
    }
}

for(i in 1:length(argweaver_trees_list)){
    avg_trajs_neut_aw[,i] <- rowMeans(trajs_neut_aw[,i,], na.rm = TRUE)
    avg_vars_neut_bin_aw[,i] <- rowMeans(vars_neut_bin_aw[,i,], na.rm = TRUE)
}

avg_trajs_neut_aw[time == 0,] <- (n_ders/n_chroms)#in the present, just use sample allele frequency.
#This is the same as the output of est_af_traj_neut() if no coalescent times get rounded to 0.
avg_vars_neut_bin_aw[time == 0,] <- (n_ders/n_chroms)*(1 - (n_ders/n_chroms))/n_chroms

traj.phen.neut.aw <- 2 * avg_trajs_neut_aw %*%  eff_sizes 
var.phen.neut.bin.aw <- 4 * avg_vars_neut_bin_aw %*% eff_sizes^2

#Method of moments from smoothed coalescent time estimates.
trajs_mom_smoothtime_aw <- array(dim = c(length(time), length(argweaver_trees_list), length(argweaver_trees_list[[1]])))# lineages remaining estimator
trajs_var_mom_smoothtime_aw <- array(dim = c(length(time), length(argweaver_trees_list), length(argweaver_trees_list[[1]])))

avg_trajs_mom_smoothtime_aw <-  matrix(nrow = length(time), ncol = length(argweaver_trees_list))
avg_trajs_var_mom_smoothtime_aw <- matrix(nrow = length(time), ncol = length(argweaver_trees_list))

for(i in 1:length(argweaver_trees_list)){	
    for(j in 21:length(argweaver_trees_list[[i]])){
    	trajs_mom_smoothtime_aw[,i,j] <- est_af_traj_mom.smoothtime(i, lins.list.argweaver[[i]][[j]], time)
	    trajs_var_mom_smoothtime_aw[,i,j] <- est_af_var_mom.smoothtime(lins.list.argweaver[[i]][[j]], time*2*N)
    }
}

for(i in 1:length(argweaver_trees_list)){
    avg_trajs_mom_smoothtime_aw[,i] <- rowMeans(trajs_mom_smoothtime_aw[,i,], na.rm = TRUE)
    avg_trajs_var_mom_smoothtime_aw[,i] <- rowMeans(trajs_var_mom_smoothtime_aw[,i,], na.rm = TRUE)
}


traj.phen.mom_smoothtime_aw <- 2 * avg_trajs_mom_smoothtime_aw %*%  eff_sizes 
var.phen.mom_smoothtime_aw <- 4 * avg_trajs_var_mom_smoothtime_aw %*%  eff_sizes^2 

traj.phen.mom_smoothtime_aw[time == 0] <- traj.phen.neut.aw[time == 0]
var.phen.mom_smoothtime_aw[time == 0] <- var.phen.neut.bin.aw[time == 0]



#waiting time-based estimates and variance
trajs_est_wt_l1_aw <- array(dim = c(length(time), length(argweaver_trees_list), length(argweaver_trees_list[[1]])))
trajs_var_wt_l1_aw <- array(dim = c(length(time), length(argweaver_trees_list), length(argweaver_trees_list[[1]])))

avg_trajs_est_wt_l1_aw <- matrix(nrow = length(time), ncol = length(argweaver_trees_list))
avg_trajs_var_wt_l1_aw <- matrix(nrow = length(time), ncol = length(argweaver_trees_list))

for(i in 1:length(argweaver_trees_list)){
    for(j in 21:length(argweaver_trees_list[[i]])){
        wt.estvar.aw <- p_ests_wait(i, anc_trees_argweaver[[i]][[j]], der_trees_argweaver[[i]][[j]], lins.list.argweaver[[i]][[j]], times.c.argweaver[[i]][[j]], time, ell.ref = 1, ell.alt = 1)
        trajs_est_wt_l1_aw[,i,j] <- wt.estvar.aw[,1]
        trajs_var_wt_l1_aw[,i,j] <- wt.estvar.aw[,2]
    }
}

for(i in 1:length(argweaver_trees_list)){
    avg_trajs_est_wt_l1_aw[,i] <- rowMeans(trajs_est_wt_l1_aw[,i,], na.rm = TRUE)
    avg_trajs_var_wt_l1_aw[,i] <- rowMeans(trajs_var_wt_l1_aw[,i,], na.rm = TRUE)
}

traj.phen.wt_l1_aw <- 2 * avg_trajs_est_wt_l1_aw %*%  eff_sizes 
var.phen.wt_l1_aw <- 4 * avg_trajs_var_wt_l1_aw %*%  eff_sizes^2 

traj.phen.wt_l1_aw[time == 0] <- traj.phen.neut.aw[time == 0]
var.phen.wt_l1_aw[time == 0] <- var.phen.neut.bin.aw[time == 0]



true.per.time <- numeric(0)
for(j in 1:length(time)){
	true.per.time[j] <- phen.traj[which.min(abs(pt.time - time[j]))]
}

mat.trajs <- matrixify.list.of.trajs(trajs)

true.afs.per.time <- matrix(nrow = length(time), ncol = n.loci)
for(j in 1:length(time)){
	true.afs.per.time[j,] <- mat.trajs[which.min(abs(pt.time - time[j])),]
}


#save errors, unscaled and scaled by estimated se, for each method at all times

err.neut.aw <- traj.phen.neut.aw - true.per.time
err.smoothmom.aw <- traj.phen.mom_smoothtime_aw - true.per.time
err.wt_l1.aw <- traj.phen.wt_l1_aw - true.per.time

err.neut.std.bin.aw <- err.neut.aw / sqrt(var.phen.neut.bin.aw)
err.smoothmom.std.aw <- err.smoothmom.aw / sqrt(var.phen.mom_smoothtime_aw)
err.wt_l1.std.aw <- err.wt_l1.aw / sqrt(var.phen.wt_l1_aw)

err.mat.aw <- cbind(err.neut.aw, err.smoothmom.aw, err.wt_l1.aw)
err.mat.std.aw <- cbind(err.neut.std.bin.aw, err.smoothmom.std.aw, err.wt_l1.std.aw)

err.array.aw[,,new_iter] <- err.mat.aw
err.std.array.aw[,,new_iter] <- err.mat.std.aw
mat.true.phentrajs[,new_iter] <- true.per.time


#qxtest_mat[new_iter,1:3] <- Qx_test(avg_trajs_neut_aw[time %in% ((0:10)/100),], eff_sizes, perms = 0)
#qxtest_mat[new_iter,4] <- Qx_test(avg_trajs_neut_aw[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
#qxtest_mat[new_iter,5:7] <- Qx_test(avg_trajs_mom_smoothtime_aw[time %in% ((0:10)/100),], eff_sizes, perms = 0)
#qxtest_mat[new_iter,8] <- Qx_test(avg_trajs_mom_smoothtime_aw[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
#qxtest_mat[new_iter,9:11] <- Qx_test(avg_trajs_est_wt_l1_aw[time %in% ((0:10)/100),], eff_sizes, perms = 0)
#qxtest_mat[new_iter,12] <- Qx_test(avg_trajs_est_wt_l1_aw[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
