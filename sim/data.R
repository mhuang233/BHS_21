# data

library(gtools)
covariates <- c("x1", "x2", "x3", "x4")
cov1 <- combinations(4, 1, covariates, repeats.allowed = F) # 4
cov2 <- combinations(4, 2, covariates, repeats.allowed = F) # 6
cov3 <- combinations(4, 3, covariates, repeats.allowed = F) # 4
cov4 <- combinations(4, 4, covariates, repeats.allowed = F) # 1

# model strings
model_strings <- NULL
var_names <- NULL

# for cov 1

for (k in 1:4){
  covx <- paste0(cov1[k, ], collapse = " + ")
  model_strings[k] <- paste0("y ~ ", covx, " + (1|i)")
  var_names[k] <- paste0(cov1[k, ], collapse = "_")
}

# for cov 2

for (k in 1:6){
  covx <- paste0(cov2[k, ], collapse = " + ")
  model_strings[k+4] <- paste0("y ~ ", covx, " + (1|i)")
  var_names[k+4] <- paste0(cov2[k, ], collapse = "_")
}

# for cov 3
for (k in 1:4){
  
  covx <- paste0(cov3[k, ], collapse = " + ")
  model_strings[k+10] <- paste0("y ~ ", covx, " + (1|i)")
  var_names[k+10] <- paste0(cov3[k, ], collapse = "_")
  
}


# for cov 4

covx <- paste0(cov4, collapse = " + ")
model_strings[15] <- paste0("y ~ ", covx, " + (1|i)")
var_names[15] <- paste0(cov4, collapse = "_")

# check model strings
model_strings
var_names

# arguements
# ni*nj: 150*30, 50*30
# icc: 0.1, 0.2, 0.3
# sigma: 0, 1, 5
# rep: 1-100

arg <- data.frame(
  ni = c(rep(150, 9), rep(50, 9)),
  nj = c(rep(30, 9), rep(10, 9)),
  icc = rep(c(0.1, 0.2, 0.3),6),
  sigma = rep(c(0, 1, 5), each = 3),
  gb = c(rep(45, 9), rep(10, 9))
)

args <- do.call("rbind", replicate(100, arg, simplify = FALSE))
reps <- rep(1:100, each = 18) 
seed <- sample(0:999, 100)
seeds <- rep(seed, each = 18)
arguement <- cbind(reps, seeds, args)

write.table(arguement, file = "arg.txt",
            col.names = F, row.names = F)

# large
arg <- read.table("arg.txt")

seed <- arg[,2] %>% unique()

# saveRDS(seed, file = "seed.rds")
argl <- data.frame(
  ni = rep(150, 9),
  nj = rep(30, 9),
  icc = rep(c(0.1, 0.2, 0.3), 3),
  sigma = rep(c(0, 5, 10), each = 3),
  gb = rep(40, 9)
)

args <- do.call("rbind", replicate(100, argl, simplify = FALSE))
reps <- rep(1:100, each = 9) 
seeds <- rep(seed, each = 9)
arguement <- cbind(reps, seeds, args) 


write.table(arguement, file = "argl.txt",
            col.names = F, row.names = F)





