#once a trait dataset is loaded, analyze all the argneedle trees.


#neutral - argneedle
trajs_neut_argneedle <- matrix(nrow = length(time), ncol = length(argneedle_trees_list))
vars_neut_bin_argneedle <- matrix(nrow = length(time), ncol = length(argneedle_trees_list))
#vars_neut_post_argneedle <- matrix(nrow = length(time), ncol = length(argneedle_trees_list))
for(i in 1:length(argneedle_trees_list)){	
  trajs_neut_argneedle[,i] <- est_af_traj_neut(lins.list.argneedle[[i]])
  vars_neut_bin_argneedle[,i] <- est_af_var_neut_bin(lins.list.argneedle[[i]])
  #	vars_neut_post_argneedle[,i] <- est_af_var_neut_post(lins.list.argneedle[[i]])
}
trajs_neut_argneedle[time == 0,] <- (n_ders/n_chroms) #in the present, just use sample allele frequency.
#This is the same as the output of est_af_traj_neut() if no coalescent times get rounded to 0.
vars_neut_bin_argneedle[time == 0,] <- (n_ders/n_chroms)*(1 - (n_ders/n_chroms))/n_chroms
traj.phen.neut.argneedle <- 2 * trajs_neut_argneedle %*%  eff_sizes 
var.phen.neut.bin.argneedle <- 4 * vars_neut_bin_argneedle %*% eff_sizes^2
#var.phen.neut.post.argneedle <- 4 * vars_neut_post_argneedle %*% eff_sizes^2

#Method of moments from smoothed coalescent time estimates---argneedle.
trajs_mom_smoothtime_argneedle <- matrix(nrow = length(time), ncol = length(argneedle_trees_list))
trajs_var_mom_smoothtime_argneedle <- matrix(nrow = length(time), ncol = length(argneedle_trees_list))
for(i in 1:length(argneedle_trees_list)){		
  trajs_mom_smoothtime_argneedle[,i] <- est_af_traj_mom.smoothtime(i, lins.list.argneedle[[i]], time)
  trajs_var_mom_smoothtime_argneedle[,i] <- est_af_var_mom.smoothtime(lins.list.argneedle[[i]], time*2*N)
}
traj.phen.mom_smoothtime_argneedle <- 2 * trajs_mom_smoothtime_argneedle %*%  eff_sizes 
var.phen.mom_smoothtime_argneedle <- 4 * trajs_var_mom_smoothtime_argneedle %*%  eff_sizes^2 
traj.phen.mom_smoothtime_argneedle[time == 0] <- traj.phen.neut.argneedle[time == 0]
var.phen.mom_smoothtime_argneedle[time == 0] <- var.phen.neut.bin.argneedle[time == 0]

#waiting time-based estimates and variance---argneedle
trajs_est_wt_l1_argneedle <- matrix(nrow = length(time), ncol = length(argneedle_trees_list))
trajs_var_wt_l1_argneedle <- matrix(nrow = length(time), ncol = length(argneedle_trees_list))
for(i in 1:length(argneedle_trees_list)){
  wt.estvar.argneedle <- p_ests_wait(i, anc_trees_argneedle[[i]], der_trees_argneedle[[i]], lins.list.argneedle[[i]], times.c.argneedle[[i]], time, ell.ref = 1, ell.alt = 1)
  trajs_est_wt_l1_argneedle[,i] <- wt.estvar.argneedle[,1]	
  trajs_var_wt_l1_argneedle[,i] <- wt.estvar.argneedle[,2]	
}
traj.phen.wt_l1_argneedle <- 2 * trajs_est_wt_l1_argneedle %*%  eff_sizes 
var.phen.wt_l1_argneedle <- 4 * trajs_var_wt_l1_argneedle %*%  eff_sizes^2 
traj.phen.wt_l1_argneedle[time == 0] <- traj.phen.neut.argneedle[time == 0]
var.phen.wt_l1_argneedle[time == 0] <- var.phen.neut.bin.argneedle[time == 0]

true.per.time <- numeric(0)
for(j in 1:length(time)){
  true.per.time[j] <- phen.traj[which.min(abs(pt.time - time[j]))]
}

mat.trajs <- matrixify.list.of.trajs(trajs)

true.afs.per.time <- matrix(nrow = length(time), ncol = n.loci)
for(j in 1:length(time)){
  true.afs.per.time[j,] <- mat.trajs[which.min(abs(pt.time - time[j])),]
}

err.neut.argneedle <- traj.phen.neut.argneedle - true.per.time
err.smoothmom.argneedle <- traj.phen.mom_smoothtime_argneedle - true.per.time
err.wt_l1.argneedle <- traj.phen.wt_l1_argneedle - true.per.time
#err.sharedN.wt.argneedle <- traj.phen.sharedN.wt.argneedle - true.per.time

err.neut.std.bin.argneedle <- err.neut.argneedle / sqrt(var.phen.neut.bin.argneedle)
err.smoothmom.std.argneedle <- err.smoothmom.argneedle / sqrt(var.phen.mom_smoothtime_argneedle)
err.wt_l1.std.argneedle <- err.wt_l1.argneedle / sqrt(var.phen.wt_l1_argneedle)
#err.sharedN.wt.std.argneedle <- err.sharedN.wt.argneedle / sqrt(var.phen.sharedN.wt.argneedle)

err.mat.argneedle <- cbind(err.neut.argneedle, err.smoothmom.argneedle, err.wt_l1.argneedle)
err.mat.std.argneedle <- cbind(err.neut.std.bin.argneedle, err.smoothmom.std.argneedle, err.wt_l1.std.argneedle)

err.array.an[,,new_iter] <- err.mat.argneedle
err.std.array.an[,,new_iter] <- err.mat.std.argneedle
mat.true.phentrajs[,new_iter] <- true.per.time

#Tests
qxtest_mat_argneedle[new_iter,1:3] <- Qx_test(trajs_neut_argneedle[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_argneedle[new_iter,4] <- Qx_test(trajs_neut_argneedle[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat_argneedle[new_iter,5:7] <- Qx_test(trajs_mom_smoothtime_argneedle[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_argneedle[new_iter,8] <- Qx_test(trajs_mom_smoothtime_argneedle[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat_argneedle[new_iter,9:11] <- Qx_test(trajs_est_wt_l1_argneedle[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_argneedle[new_iter,12] <- Qx_test(trajs_est_wt_l1_argneedle[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
#qxtest_mat_argneedle[iter, 13:15] <- Qx_test(trajs_est_sharedN.wt.argneedle[time %in% ((0:10)/100),], eff_sizes, perms = 0)
#qxtest_mat_argneedle[iter, 16] <- Qx_test(trajs_est_sharedN.wt.argneedle[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]


#print("argneedle tree T_X statistic, number of timepoints, and permutation p")
#print(Qx_test(trajs_neut_argneedle[time %in% ((0:10)/100),], eff_sizes, perms = 10000))

