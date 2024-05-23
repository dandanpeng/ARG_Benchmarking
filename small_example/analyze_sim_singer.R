#neutral
trajs_neut_singer <- array(dim = c(length(time), length(singer_trees_list), length(singer_trees_list[[1]]))) #matrix of proportion of derived lineages (1st estimator)
vars_neut_bin_singer <- array(dim = c(length(time), length(singer_trees_list), length(singer_trees_list[[1]])))

avg_trajs_neut_singer <- matrix(nrow = length(time), ncol = length(singer_trees_list))
avg_vars_neut_bin_singer <- matrix(nrow = length(time), ncol = length(singer_trees_list))

for(i in 1:length(singer_trees_list)){	
  for(j in 1:length(singer_trees_list[[i]])){
    trajs_neut_singer[,i,j] <- est_af_traj_neut(lins.list.singer[[i]][[j]])
    vars_neut_bin_singer[,i,j] <- est_af_var_neut_bin(lins.list.singer[[i]][[j]])
  }
}

for(i in 1:length(singer_trees_list)){
  avg_trajs_neut_singer[,i] <- rowMeans(trajs_neut_singer[,i,])
  avg_vars_neut_bin_singer[,i] <- rowMeans(vars_neut_bin_singer[,i,])
}

avg_trajs_neut_singer[time == 0,] <- (n_ders/n_chroms)#in the present, just use sample allele frequency.
#This is the same as the output of est_af_traj_neut() if no coalescent times get rounded to 0.
avg_vars_neut_bin_singer[time == 0,] <- (n_ders/n_chroms)*(1 - (n_ders/n_chroms))/n_chroms

traj.phen.neut.singer <- 2 * avg_trajs_neut_singer %*%  eff_sizes 
var.phen.neut.bin.singer <- 4 * avg_vars_neut_bin_singer %*% eff_sizes^2

#Method of moments from smoothed coalescent time estimates.
trajs_mom_smoothtime_singer <- array(dim = c(length(time), length(singer_trees_list), length(singer_trees_list[[1]])))# lineages remaining estimator
trajs_var_mom_smoothtime_singer <- array(dim = c(length(time), length(singer_trees_list), length(singer_trees_list[[1]])))

avg_trajs_mom_smoothtime_singer <-  matrix(nrow = length(time), ncol = length(singer_trees_list))
avg_trajs_var_mom_smoothtime_singer <- matrix(nrow = length(time), ncol = length(singer_trees_list))

for(i in 1:length(singer_trees_list)){	
  for(j in 1:length(singer_trees_list[[i]])){
    trajs_mom_smoothtime_singer[,i,j] <- est_af_traj_mom.smoothtime(i, lins.list.singer[[i]][[j]], time)
    trajs_var_mom_smoothtime_singer[,i,j] <- est_af_var_mom.smoothtime(lins.list.singer[[i]][[j]], time*2*N)
  }
}

for(i in 1:length(singer_trees_list)){
  avg_trajs_mom_smoothtime_singer[,i] <- rowMeans(trajs_mom_smoothtime_singer[,i,])
  avg_trajs_var_mom_smoothtime_singer[,i] <- rowMeans(trajs_var_mom_smoothtime_singer[,i,])
}


traj.phen.mom_smoothtime_singer <- 2 * avg_trajs_mom_smoothtime_singer %*%  eff_sizes 
var.phen.mom_smoothtime_singer <- 4 * avg_trajs_var_mom_smoothtime_singer %*%  eff_sizes^2 

traj.phen.mom_smoothtime_singer[time == 0] <- traj.phen.neut.singer[time == 0]
var.phen.mom_smoothtime_singer[time == 0] <- var.phen.neut.bin.singer[time == 0]



#waiting time-based estimates and variance
trajs_est_wt_l1_singer <- array(dim = c(length(time), length(singer_trees_list), length(singer_trees_list[[1]])))
trajs_var_wt_l1_singer <- array(dim = c(length(time), length(singer_trees_list), length(singer_trees_list[[1]])))

avg_trajs_est_wt_l1_singer <- matrix(nrow = length(time), ncol = length(singer_trees_list))
avg_trajs_var_wt_l1_singer <- matrix(nrow = length(time), ncol = length(singer_trees_list))

for(i in 1:length(singer_trees_list)){
  print(i)
  for(j in 1:length(singer_trees_list[[i]])){
    wt.estvar.singer <- p_ests_wait(i, anc_trees_singer[[i]][[j]], der_trees_singer[[i]][[j]], lins.list.singer[[i]][[j]], times.c.singer[[i]][[j]], time, ell.ref = 1, ell.alt = 1)
    trajs_est_wt_l1_singer[,i,j] <- wt.estvar.singer[,1]
    trajs_var_wt_l1_singer[,i,j] <- wt.estvar.singer[,2]
  }
}

for(i in 1:length(singer_trees_list)){
  avg_trajs_est_wt_l1_singer[,i] <- rowMeans(trajs_est_wt_l1_singer[,i,])
  avg_trajs_var_wt_l1_singer[,i] <- rowMeans(trajs_var_wt_l1_singer[,i,])
}

traj.phen.wt_l1_singer <- 2 * avg_trajs_est_wt_l1_singer %*%  eff_sizes 
var.phen.wt_l1_singer <- 4 * avg_trajs_var_wt_l1_singer %*%  eff_sizes^2 

traj.phen.wt_l1_singer[time == 0] <- traj.phen.neut.singer[time == 0]
var.phen.wt_l1_singer[time == 0] <- var.phen.neut.bin.singer[time == 0]



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

err.neut.singer <- traj.phen.neut.singer - true.per.time
err.smoothmom.singer <- traj.phen.mom_smoothtime_singer - true.per.time
err.wt_l1.singer <- traj.phen.wt_l1_singer - true.per.time

err.neut.std.bin.singer <- err.neut.singer / sqrt(var.phen.neut.bin.singer)
err.smoothmom.std.singer <- err.smoothmom.singer / sqrt(var.phen.mom_smoothtime_singer)
err.wt_l1.std.singer <- err.wt_l1.singer / sqrt(var.phen.wt_l1_singer)

err.mat.singer <- cbind(err.neut.singer, err.smoothmom.singer, err.wt_l1.singer)
err.mat.std.singer <- cbind(err.neut.std.bin.singer, err.smoothmom.std.singer, err.wt_l1.std.singer)

err.array.singer[,,new_iter] <- err.mat.singer
err.std.array.singer[,,new_iter] <- err.mat.std.singer
mat.true.phentrajs[,new_iter] <- true.per.time


qxtest_mat_singer[new_iter,1:3] <- Qx_test(avg_trajs_neut_singer[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_singer[new_iter,4] <- Qx_test(avg_trajs_neut_singer[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat_singer[new_iter,5:7] <- Qx_test(avg_trajs_mom_smoothtime_singer[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_singer[new_iter,8] <- Qx_test(avg_trajs_mom_smoothtime_singer[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
qxtest_mat_singer[new_iter,9:11] <- Qx_test(avg_trajs_est_wt_l1_singer[time %in% ((0:10)/100),], eff_sizes, perms = 0)
qxtest_mat_singer[new_iter,12] <- Qx_test(avg_trajs_est_wt_l1_singer[time %in% ((0:10)/100),], eff_sizes, perms = 10000)[3]
