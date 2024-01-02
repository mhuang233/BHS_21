# arguement
library(tidyverse)
seed <- sample(1:1000, 100)

# for sigma = 0.1
sd1 <- sqrt(.001/.9) %>% round(3)
sd2 <- sqrt(.001/.6) %>% round(3)
sd3 <- sqrt(.001/.1) %>% round(3)

# for sigma = 5
sd4 <- sqrt(2.5/.9) %>% round(3)
sd5 <- sqrt(2.5/.6) %>% round(3)
sd6 <- sqrt(2.5/.1) %>% round(3)

# for sigma = 10
sd7 <- sqrt(10/.9) %>% round(3)
sd8 <- sqrt(10/.6) %>% round(3)
sd9 <- sqrt(10/.1) %>% round(3)

arg <- data.frame(rep = rep(1:100, each = 18),
                  c = rep(seed, each = 18),
                  ni = rep(c(50, 150), each = 9),
                  nj = rep(c(10, 30), each = 9),
                  sd = rep(c(sd4, sd5, sd6, sd4, sd5, sd6, sd7, sd8, sd9), 2), 
                  sigma = rep(c(0, 5, 10), each = 3), 
                  gb = rep(c(30, 40), each = 9))

write.table(arg, "arg_bhs.txt", col.names = F, row.names = F)


arg <- data.frame(rep = rep(1:100, each = 12),
                  c = rep(seed, each = 12),
                  ni = rep(c(50, 150), each = 6),
                  nj = rep(c(10, 30), each = 6),
                  sd = rep(c(sd4, sd5, sd6), 4), 
                  sigma = rep(c(0, 5), each = 3), 
                  gb = rep(c(35, 45), each = 6))

