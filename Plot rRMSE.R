library(ggplot2)
library(gridExtra)

RMSE_2.0 <- read.csv("RMSE2.0.csv", header=TRUE)
RMSE_2.0_partial <- read.csv("RMSE2.0_partial.csv", header=TRUE)

p1 <- ggplot(RMSE_2.0, aes(x = rRMSE, y = CI, shape = Type, color = Imputation)) +
  geom_point(size = 2) +  
  theme_minimal() +  
  labs(x = "rRMSE",
       y = "Coverage",
       color = "Imputation",
       shape = "Type") +
  scale_color_manual(values = c("blue", "red", "green", "purple", "orange")) + 
  theme(legend.title = element_text(size = 12)) + 
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) 

p2 <- ggplot(RMSE_2.0_partial, aes(x = rRMSE, y = CI, shape = Type, color = Imputation)) +
  geom_point(size = 2) +  
  theme_minimal() +  
  labs(x = "rRMSE",
       y = "Coverage",
       color = "Imputation",
       shape = "Type") +
  scale_color_manual(values = c("green", "purple", "orange", "blue", "red")) + 
  theme(legend.title = element_text(size = 12)) + 
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) 


options(repr.plot.width = 10, repr.plot.height = 8, repr.plot.res = 200)
# p1 with axis lines only for x >= 0 and y >= 0
p1 <- p1 +
  geom_segment(aes(x = 0, xend = max(RMSE_2.0$rRMSE), y = 0, yend = 0),
               inherit.aes = FALSE, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(RMSE_2.0$CI)),
               inherit.aes = FALSE, color = "black") +
  theme(panel.grid = element_blank())

# p2 with axis lines only for x >= 0 and y >= 0
p2 <- p2 +
  geom_segment(aes(x = 0, xend = max(RMSE_2.0_partial$rRMSE), y = 80, yend = 80),
               inherit.aes = FALSE, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 80, yend = max(RMSE_2.0_partial$CI)),
               inherit.aes = FALSE, color = "black") +
  theme(panel.grid = element_blank())



