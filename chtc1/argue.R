# arguement
#arg1 <- read.table("argl.txt")
#c <- arg1$V2 %>% unique()

seed <- sample(1:1000, 100)

# for icc = 1
sd1 <- sqrt(2.5/.9) %>% round(3)
sd2 <- sqrt(2.5/.6) %>% round(3)
sd3 <- sqrt(2.5/.1) %>% round(3)


arg <- data.frame(rep = rep(1:100, each = 18),
                  c = rep(seed, each = 18),
                  ni = rep(c(50, 150), each = 9),
                  nj = rep(c(10, 30), each = 9),
                  sd = rep(c(sd1, sd2, 5.000), 6), 
                  sigma = rep(c(0, 5, 10), each = 3), 
                  gb = rep(c(30, 35, 40, 40, 45, 55), each = 3))

write.table(arg, "arg_bhs.txt", col.names = F, row.names = F)
