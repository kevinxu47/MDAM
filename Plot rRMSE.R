library(ggplot2)


RMSE_0.5 <- read.csv("RMSE0.5.csv", header=TRUE)

ggplot(RMSE_0.5, aes(x = rRMSE, y = CI, shape = Type, color = Imputation)) +
  geom_point(size = 1) +  # Adjust size of points
  theme_minimal() + 
  labs(title = "Coverage vs rRMSE",
       x = "rRMSE",
       y = "Coverage",
       color = "Imputation",
       shape = "Type")  +
  scale_color_manual(values = c("blue", "red", "green", "purple", "orange")) + 
  theme(legend.title = element_text(size = 12)) 




RMSE_2.0 <- read.csv("RMSE2.0.csv", header=TRUE)

ggplot(RMSE_2.0, aes(x = rRMSE, y = CI, shape = Type, color = Imputation)) +
  geom_point(size = 1) +  
  theme_minimal() +  
  labs(title = "Coverage vs rRMSE",
       x = "rRMSE",
       y = "Coverage",
       color = "Imputation",
       shape = "Type") +
  scale_color_manual(values = c("blue", "red", "green", "purple", "orange")) + 
  theme(legend.title = element_text(size = 12))  




RMSE_0.5_partial <- read.csv("RMSE0.5_partial.csv", header=TRUE)

ggplot(RMSE_0.5_partial, aes(x = rRMSE, y = CI, shape = Type, color = Imputation)) +
  geom_point(size = 1) +  
  theme_minimal() +  
  labs(title = "Coverage vs rRMSE",
       x = "rRMSE",
       y = "Coverage",
       color = "Imputation",
       shape = "Type") +
  scale_color_manual(values = c("green", "purple", "orange", "blue", "red")) + 
  theme(legend.title = element_text(size = 12)) 




RMSE_2.0_partial <- read.csv("RMSE2.0_partial.csv", header=TRUE)

ggplot(RMSE_2.0_partial, aes(x = rRMSE, y = CI, shape = Type, color = Imputation)) +
  geom_point(size = 1) +  
  theme_minimal() +  
  labs(title = "Coverage vs rRMSE",
       x = "rRMSE",
       y = "Coverage",
       color = "Imputation",
       shape = "Type") +
  scale_color_manual(values = c("green", "purple", "orange", "blue", "red")) + 
  theme(legend.title = element_text(size = 12)) 




