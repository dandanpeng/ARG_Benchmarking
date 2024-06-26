#once a trait dataset is loaded, analyze all the relate trees.


#neutral - relate
trajs_neut_relate <- matrix(nrow = length(time), ncol = length(relate_trees_list))
vars_neut_bin_relate <- matrix(nrow = length(time), ncol = length(relate_trees_list))
#vars_neut_post_relate <- matrix(nrow = length(time), ncol = length(relate_trees_list))
for(i in 1:length(relate_trees_list)){	
  trajs_neut_relate[,i] <- est_af_traj_neut(lins.list.relate[[i]])
  vars_neut_bin_relate[,i] <- est_af_var_neut_bin(lins.list.relate[[i]])
  #	vars_neut_post_relate[,i] <- est_af_var_neut_post(lins.list.relate[[i]])
}
trajs_neut_relate[time == 0,] <- (n_ders/n_chroms) #in the present, just use sample allele frequency.
#This is the same as the output of est_af_traj_neut() if no coalescent times get rounded to 0.
vars_neut_bin_relate[time == 0,] <- (n_ders/n_chroms)*(1 - (n_ders/n_chroms))/n_chroms
traj.phen.neut.relate <- 2 * trajs_neut_relate %*%  eff_sizes 
var.phen.neut.bin.relate <- 4 * vars_neut_bin_relate %*% eff_sizes^2
#var.phen.neut.post.relate <- 4 * vars_neut_post_relate %*% eff_sizes^2


#Method of moments from smoothed coalescent time estimates---relate.
trajs_mom_smoothtime_relate <- matrix(nrow = length(time), ncol = length(relate_trees_list))
trajs_var_mom_smoothtime_relate <- matrix(nrow = length(time), ncol = length(relate_trees_list))
for(i in 1:length(relate_trees_list)){	
  trajs_mom_smoothtime_relate[,i] <- est_af_traj_mom.smoothtime(i, lins.list.relate[[i]], time)
  trajs_var_mom_smoothtime_relate[,i] <- est_af_var_mom.smoothtime(lins.list.relate[[i]], time*2*N)
}
traj.phen.mom_smoothtime_relate <- 2 * trajs_mom_smoothtime_relate %*%  eff_sizes 
var.phen.mom_smoothtime_relate <- 4 * trajs_var_mom_smoothtime_relate %*%  eff_sizes^2 
traj.phen.mom_smoothtime_relate[time == 0] <- traj.phen.neut.relate[time == 0]
var.phen.mom_smoothtime_relate[time == 0] <- var.phen.neut.bin.relate[time == 0]


#waiting time-based estimates and variance---relate
trajs_est_wt_l1_relate <- matrix(nrow = length(time), ncol = length(relate_trees_list))
trajs_var_wt_l1_relate <- matrix(nrow = length(time), ncol = length(relate_trees_list))
for(i in 1:length(relate_trees_list)){
  wt.estvar.relate <- p_ests_wait(i, anc_trees_relate[[i]], der_trees_relate[[i]], lins.list.relate[[i]], times.c.relate[[i]], time, ell.ref = 1, ell.alt = 1)
  trajs_est_wt_l1_relate[,i] <- wt.estvar.relate[,1]	
  trajs_var_wt_l1_relate[,i] <- wt.estvar.relate[,2]	
}
traj.phen.wt_l1_relate <- 2 * trajs_est_wt_l1_relate %*%  eff_sizes 
var.phen.wt_l1_relate <- 4 * trajs_var_wt_l1_relate %*%  eff_sizes^2 
traj.phen.wt_l1_relate[time == 0] <- traj.phen.neut.relate[time == 0]
var.phen.wt_l1_relate[time == 0] <- var.phen.neut.bin.relate[time == 0]
print("waittime succeed!")

true.per.time <- numeric(0)
for(j in 1:length(time)){
  true.per.time[j] <- phen.traj[which.min(abs(pt.time - time[j]))]
}

mat.trajs <- matrixify.list.of.trajs(trajs)

true.afs.per.time <- matrix(nrow = length(time), ncol = n.loci)
for(j in 1:length(time)){
  true.afs.per.time[j,] <- mat.trajs[which.min(abs(pt.time - time[j])),]
}

err.neut.relate <- traj.phen.neut.relate - true.per.time
err.smoothmom.relate <- traj.phen.mom_smoothtime_relate - true.per.time
err.wt_l1.relate <- traj.phen.wt_l1_relate - true.per.time
#err.sharedN.wt.relate <- traj.phen.sharedN.wt.relate - true.per.time

err.neut.std.bin.relate <- err.neut.relate / sqrt(var.phen.neut.bin.relate)
err.smoothmom.std.relate <- err.smoothmom.relate / sqrt(var.phen.mom_smoothtime_relate)
err.wt_l1.std.relate <- err.wt_l1.relate / sqrt(var.phen.wt_l1_relate)
#err.sharedN.wt.std.relate <- err.sharedN.wt.relate / sqrt(var.phen.sharedN.wt.relate)

err.mat.relate <- cbind(err.neut.relate, err.smoothmom.relate, err.wt_l1.relate)
err.mat.std.relate <- cbind(err.neut.std.bin.relate, err.smoothmom.std.relate, err.wt_l1.std.relate)

err.array.relate[,,iteration] <- err.mat.relate
err.std.array.relate[,,iteration] <- err.mat.std.relate
mat.true.phentrajs[,iteration] <- true.per.time

#Tests
qxtest_mat_relate[iteration,1:3] <- Qx_test(trajs_neut_relate[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_relate[iteration,4] <- Qx_test(trajs_neut_relate[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat_relate[iteration,5:7] <- Qx_test(trajs_mom_smoothtime_relate[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_relate[iteration,8] <- Qx_test(trajs_mom_smoothtime_relate[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat_relate[iteration,9:11] <- Qx_test(trajs_est_wt_l1_relate[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_relate[iteration,12] <- Qx_test(trajs_est_wt_l1_relate[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
#qxtest_mat_relate[iter, 13:15] <- Qx_test(trajs_est_sharedN.wt.relate[time %in% ((0:10)/100),], eff_sizes, perms = 0)
#qxtest_mat_relate[iter, 16] <- Qx_test(trajs_est_sharedN.wt.relate[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]


print("relate tree T_X statistic, number of timepoints, and permutation p")
print(Qx_test(trajs_neut_relate[time %in% ((0:10)/100),], eff_sizes, perms = 10000))

