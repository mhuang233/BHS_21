# arguement
arg1 <- read.table("argl.txt")
c <- arg1$V2 %>% unique()

seed <- sample(1:1000, 100)

# for icc = 1
sd1 <- sqrt(2.5/.9)
sd2 <- sqrt(2.5/.6)
sd3 <- sqrt(2.5/.1)


arg <- data.frame(rep = rep(1:100, each = 12),
                  c = rep(seed, each = 12),
                  ni = rep(c(50, 150), each = 6),
                  nj = rep(c(10, 30), each = 6),
                  sd = rep(c(sd1, sd2, sd3), 4), 
                  sigma = rep(c(0, 5), each = 6), 
                  gb = rep(c(35, 45), each = 6))

write.table(arg, "arg_bhs.txt", col.names = F, row.names = F)
