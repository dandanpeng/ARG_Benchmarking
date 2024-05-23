
#Doc Edge, 6/26/18

#once a trait dataset is loaded, analyze all the tsinfer trees.


#neutral - tsinfer
trajs_neut_tsinfer <- matrix(nrow = length(time), ncol = length(tsinfer_trees_list))
vars_neut_bin_tsinfer <- matrix(nrow = length(time), ncol = length(tsinfer_trees_list))
#vars_neut_post_tsinfer <- matrix(nrow = length(time), ncol = length(tsinfer_trees_list))
for(i in 1:length(tsinfer_trees_list)){	
	trajs_neut_tsinfer[,i] <- est_af_traj_neut(lins.list.tsinfer[[i]])
	vars_neut_bin_tsinfer[,i] <- est_af_var_neut_bin(lins.list.tsinfer[[i]])
#	vars_neut_post_tsinfer[,i] <- est_af_var_neut_post(lins.list.tsinfer[[i]])
}
trajs_neut_tsinfer[time == 0,] <- (n_ders/n_chroms) #in the present, just use sample allele frequency.
#This is the same as the output of est_af_traj_neut() if no coalescent times get rounded to 0.
vars_neut_bin_tsinfer[time == 0,] <- (n_ders/n_chroms)*(1 - (n_ders/n_chroms))/n_chroms
traj.phen.neut.tsinfer <- 2 * trajs_neut_tsinfer %*%  eff_sizes 
var.phen.neut.bin.tsinfer <- 4 * vars_neut_bin_tsinfer %*% eff_sizes^2
#var.phen.neut.post.tsinfer <- 4 * vars_neut_post_tsinfer %*% eff_sizes^2




#Method of moments from smoothed coalescent time estimates---tsinfer.
trajs_mom_smoothtime_tsinfer <- matrix(nrow = length(time), ncol = length(tsinfer_trees_list))
trajs_var_mom_smoothtime_tsinfer <- matrix(nrow = length(time), ncol = length(tsinfer_trees_list))
for(i in 1:length(tsinfer_trees_list)){		
	trajs_mom_smoothtime_tsinfer[,i] <- est_af_traj_mom.smoothtime(i, lins.list.tsinfer[[i]], time)
	trajs_var_mom_smoothtime_tsinfer[,i] <- est_af_var_mom.smoothtime(lins.list.tsinfer[[i]], time*2*N)
}
traj.phen.mom_smoothtime_tsinfer <- 2 * trajs_mom_smoothtime_tsinfer %*%  eff_sizes 
var.phen.mom_smoothtime_tsinfer <- 4 * trajs_var_mom_smoothtime_tsinfer %*%  eff_sizes^2 
traj.phen.mom_smoothtime_tsinfer[time == 0] <- traj.phen.neut.tsinfer[time == 0]
var.phen.mom_smoothtime_tsinfer[time == 0] <- var.phen.neut.bin.tsinfer[time == 0]


#waiting time-based estimates and variance---tsinfer
trajs_est_wt_l1_tsinfer <- matrix(nrow = length(time), ncol = length(tsinfer_trees_list))
trajs_var_wt_l1_tsinfer <- matrix(nrow = length(time), ncol = length(tsinfer_trees_list))
for(i in 1:length(tsinfer_trees_list)){
	wt.estvar.tsinfer <- p_ests_wait(i, anc_trees_tsinfer[[i]], der_trees_tsinfer[[i]], lins.list.tsinfer[[i]], times.c.tsinfer[[i]], time, ell.ref = 1, ell.alt = 1)
	trajs_est_wt_l1_tsinfer[,i] <- wt.estvar.tsinfer[,1]	
	trajs_var_wt_l1_tsinfer[,i] <- wt.estvar.tsinfer[,2]	
}
traj.phen.wt_l1_tsinfer <- 2 * trajs_est_wt_l1_tsinfer %*%  eff_sizes 
var.phen.wt_l1_tsinfer <- 4 * trajs_var_wt_l1_tsinfer %*%  eff_sizes^2 
traj.phen.wt_l1_tsinfer[time == 0] <- traj.phen.neut.tsinfer[time == 0]
var.phen.wt_l1_tsinfer[time == 0] <- var.phen.neut.bin.tsinfer[time == 0]



true.per.time <- numeric(0)
for(j in 1:length(time)){
	true.per.time[j] <- phen.traj[which.min(abs(pt.time - time[j]))]
}

mat.trajs <- matrixify.list.of.trajs(trajs)

true.afs.per.time <- matrix(nrow = length(time), ncol = n.loci)
for(j in 1:length(time)){
	true.afs.per.time[j,] <- mat.trajs[which.min(abs(pt.time - time[j])),]
}



#save errors, unscaled and scaled by estimated se, for each method at all times---tsinfer.

err.neut.tsinfer <- traj.phen.neut.tsinfer - true.per.time
err.smoothmom.tsinfer <- traj.phen.mom_smoothtime_tsinfer - true.per.time
err.wt_l1.tsinfer <- traj.phen.wt_l1_tsinfer - true.per.time


err.neut.std.bin.tsinfer <- err.neut.tsinfer / sqrt(var.phen.neut.bin.tsinfer)
err.smoothmom.std.tsinfer <- err.smoothmom.tsinfer / sqrt(var.phen.mom_smoothtime_tsinfer)
err.wt_l1.std.tsinfer <- err.wt_l1.tsinfer / sqrt(var.phen.wt_l1_tsinfer)

err.mat.tsinfer <- cbind(err.neut.tsinfer, err.smoothmom.tsinfer, err.wt_l1.tsinfer)
err.mat.std.tsinfer <- cbind(err.neut.std.bin.tsinfer, err.smoothmom.std.tsinfer, err.wt_l1.std.tsinfer)

err.array.tsinfer[,,iteration] <- err.mat.tsinfer
err.std.array.tsinfer[,,iteration] <- err.mat.std.tsinfer
mat.true.phentrajs[,iteration] <- true.per.time


#Tests
qxtest_mat_tsinfer[iteration,1:3] <- Qx_test(trajs_neut_tsinfer[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_tsinfer[iteration,4] <- Qx_test(trajs_neut_tsinfer[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat_tsinfer[iteration,5:7] <- Qx_test(trajs_mom_smoothtime_tsinfer[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_tsinfer[iteration,8] <- Qx_test(trajs_mom_smoothtime_tsinfer[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat_tsinfer[iteration,9:11] <- Qx_test(trajs_est_wt_l1_tsinfer[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_tsinfer[iteration,12] <- Qx_test(trajs_est_wt_l1_tsinfer[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]

#print("tsinfer tree T_X statistic, number of timepoints, and permutation p")
#print(Qx_test(trajs_neut_tsinfer[time %in% ((0:10)/100),], eff_sizes, perms = 10000))

