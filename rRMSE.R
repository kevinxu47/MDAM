
# Calculate rRMSE for cond_prob
rRMSE_cond_prob_1 <- (1 / cond_prob_true[1]) * sqrt(sum((cond_prob[,1] - cond_prob_true[1])^2) / (nrow(cond_prob) - 1))
rRMSE_cond_prob_2 <- (1 / cond_prob_true[2]) * sqrt(sum((cond_prob[,2] - cond_prob_true[2])^2) / (nrow(cond_prob) - 1))
rRMSE_cond_prob_3 <- (1 / cond_prob_true[3]) * sqrt(sum((cond_prob[,3] - cond_prob_true[3])^2) / (nrow(cond_prob) - 1))
rRMSE_cond_prob_4 <- (1 / cond_prob_true[4]) * sqrt(sum((cond_prob[,4] - cond_prob_true[4])^2) / (nrow(cond_prob) - 1))
rRMSE_cond_prob_5 <- (1 / cond_prob_true[5]) * sqrt(sum((cond_prob[,5] - cond_prob_true[5])^2) / (nrow(cond_prob) - 1))
rRMSE_cond_prob_6 <- (1 / cond_prob_true[6]) * sqrt(sum((cond_prob[,6] - cond_prob_true[6])^2) / (nrow(cond_prob) - 1))
rRMSE_cond_prob_7 <- (1 / cond_prob_true[7]) * sqrt(sum((cond_prob[,7] - cond_prob_true[7])^2) / (nrow(cond_prob) - 1))
rRMSE_cond_prob_8 <- (1 / cond_prob_true[8]) * sqrt(sum((cond_prob[,8] - cond_prob_true[8])^2) / (nrow(cond_prob) - 1))

# Calculate rRMSE for joint_prob
rRMSE_joint_prob_1 <- (1 / joint_prob_true[1]) * sqrt(sum((joint_prob[,1] - joint_prob_true[1])^2) / (nrow(joint_prob) - 1))
rRMSE_joint_prob_2 <- (1 / joint_prob_true[2]) * sqrt(sum((joint_prob[,2] - joint_prob_true[2])^2) / (nrow(joint_prob) - 1))
rRMSE_joint_prob_3 <- (1 / joint_prob_true[3]) * sqrt(sum((joint_prob[,3] - joint_prob_true[3])^2) / (nrow(joint_prob) - 1))
rRMSE_joint_prob_4 <- (1 / joint_prob_true[4]) * sqrt(sum((joint_prob[,4] - joint_prob_true[4])^2) / (nrow(joint_prob) - 1))

# Calculate rRMSE for cond_prob_given12
rRMSE_cond_prob_given12_1 <- (1 / cond_prob_given12_true[1]) * sqrt(sum((cond_prob_given12[,1] - cond_prob_given12_true[1])^2) / (nrow(cond_prob_given12) - 1))
rRMSE_cond_prob_given12_2 <- (1 / cond_prob_given12_true[2]) * sqrt(sum((cond_prob_given12[,2] - cond_prob_given12_true[2])^2) / (nrow(cond_prob_given12) - 1))
rRMSE_cond_prob_given12_3 <- (1 / cond_prob_given12_true[3]) * sqrt(sum((cond_prob_given12[,3] - cond_prob_given12_true[3])^2) / (nrow(cond_prob_given12) - 1))
rRMSE_cond_prob_given12_4 <- (1 / cond_prob_given12_true[4]) * sqrt(sum((cond_prob_given12[,4] - cond_prob_given12_true[4])^2) / (nrow(cond_prob_given12) - 1))
rRMSE_cond_prob_given12_5 <- (1 / cond_prob_given12_true[5]) * sqrt(sum((cond_prob_given12[,5] - cond_prob_given12_true[5])^2) / (nrow(cond_prob_given12) - 1))
rRMSE_cond_prob_given12_6 <- (1 / cond_prob_given12_true[6]) * sqrt(sum((cond_prob_given12[,6] - cond_prob_given12_true[6])^2) / (nrow(cond_prob_given12) - 1))
rRMSE_cond_prob_given12_7 <- (1 / cond_prob_given12_true[7]) * sqrt(sum((cond_prob_given12[,7] - cond_prob_given12_true[7])^2) / (nrow(cond_prob_given12) - 1))
rRMSE_cond_prob_given12_8 <- (1 / cond_prob_given12_true[8]) * sqrt(sum((cond_prob_given12[,8] - cond_prob_given12_true[8])^2) / (nrow(cond_prob_given12) - 1))

# Calculate rRMSE for T_est
T_X1_true <- sum(population$X1)
T_X2_true <- sum(population$X2)
T_Y1_true <- sum(population$Y1)
T_Y2_true <- sum(population$Y2)
T_Y3_true <- sum(population$Y3)
T_Y4_true <- sum(population$Y4)

rRMSE_T_est_1 <- (1 / T_X1_true) * sqrt(sum((T_est[,1] - T_X1_true)^2) / (nrow(T_est) - 1))
rRMSE_T_est_2 <- (1 / T_X2_true) * sqrt(sum((T_est[,2] - T_X2_true)^2) / (nrow(T_est) - 1))
rRMSE_T_est_3 <- (1 / T_Y1_true) * sqrt(sum((T_est[,3] - T_Y1_true)^2) / (nrow(T_est) - 1))
rRMSE_T_est_4 <- (1 / T_Y2_true) * sqrt(sum((T_est[,4] - T_Y2_true)^2) / (nrow(T_est) - 1))
rRMSE_T_est_5 <- (1 / T_Y3_true) * sqrt(sum((T_est[,5] - T_Y3_true)^2) / (nrow(T_est) - 1))
rRMSE_T_est_6 <- (1 / T_Y4_true) * sqrt(sum((T_est[,6] - T_Y4_true)^2) / (nrow(T_est) - 1))






