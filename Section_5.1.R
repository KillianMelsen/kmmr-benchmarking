# Libraries:
suppressMessages(library(sommer))
suppressMessages(library(asreml))
suppressMessages(library(kmmr))
suppressMessages(library(Matrix))
suppressMessages(library(microbenchmark))
asreml.options(workspace = "5000mb")

# Loading data:
d <- readRDS("data/DROPS_weights.rds")
Ki.asr <- solve(readRDS("data/DROPS_K.rds")[levels(d$Gen), levels(d$Gen)])
Ki.smr <- Ki.kmr <- Matrix(Ki.asr)
attr(Ki.asr, "INVERSE") <- TRUE
attr(Ki.smr, "inverse") <- TRUE
d$wt.kmr <- 1 / d$weight

# FA1/FA3 fits ====
Sys.getenv("MKL_NUM_THREADS")
## FA1 asreml ====
bm.asr <- microbenchmark(fit.asr <- asreml(fixed = GY ~ -1 + ManEnv,
                                           random = ~ fa(ManEnv, 1):vm(Gen, Ki.asr),
                                           residual = ~ units,
                                           weights = weight, family = asr_gaussian(dispersion = 1.0),
                                           data = d, maxit = 150),
                         times = 5)
tolLL <- fit.asr$trace[1, ncol(fit.asr$trace) - 1] - fit.asr$trace[1, ncol(fit.asr$trace) - 2]

## FA1 sommer ====
W <- Matrix(diag(d$weight))
bm.smr <- microbenchmark(fit.smr <- mmes(fixed = GY ~ -1 + ManEnv,
                                         random = ~ vsm(fam(ManEnv, 1), ism(Gen), Gu = Ki.smr),
                                         rcov = ~ vsm(ism(units), sigma2 = 1, fixedSigma2 = TRUE),
                                         W = W,
                                         data = d, nIters = 150, tolParConvLL = tolLL),
                         times = 5)

## FA1 kmmr ====
sv <- kmmr(fixed = GY ~ -1 + ManEnv,
                      random = ~ re(~ exfa(ManEnv, 1):lk(Gen, Ki.kmr), vars = ~ ManEnv),
                      residual = ~ units,
                      weights = "wt.kmr",
                      data = d, maxit = 150, getStartingValues = TRUE)
sv$Init[nrow(sv)] <- 1.0
sv$Fixed[nrow(sv)] <- TRUE
bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ -1 + ManEnv,
                                         random = ~ re(~ exfa(ManEnv, 1):lk(Gen, Ki.kmr), vars = ~ ManEnv),
                                         residual = ~ units,
                                         weights = "wt.kmr",
                                         data = d, maxit = 150,
                                         init = sv$Init, isFixed = sv$Fixed, tolConvLL = tolLL, denseThreshold = 0.1),
                         times = 5)

## FA3 asreml ====
bm.asr <- microbenchmark(fit.asr <- asreml(fixed = GY ~ -1 + ManEnv,
                                           random = ~ fa(ManEnv, 3):vm(Gen, Ki.asr),
                                           residual = ~ units,
                                           weights = weight, family = asr_gaussian(dispersion = 1.0),
                                           data = d, maxit = 150),
                         times = 5)
tolLL <- fit.asr$trace[1, ncol(fit.asr$trace) - 1] - fit.asr$trace[1, ncol(fit.asr$trace) - 2]

## FA3 sommer ====
W <- Matrix(diag(d$weight))
bm.smr <- microbenchmark(fit.smr <- mmes(fixed = GY ~ -1 + ManEnv,
                                         random = ~ vsm(fam(ManEnv, 3), ism(Gen), Gu = Ki.smr),
                                         rcov = ~ vsm(ism(units), sigma2 = 1, fixedSigma2 = TRUE),
                                         W = W,
                                         data = d, nIters = 150, tolParConvLL = tolLL),
                         times = 5)

## FA3 kmmr ====
sv <- kmmr(fixed = GY ~ -1 + ManEnv,
                      random = ~ re(~ exfa(ManEnv, 3):lk(Gen, Ki.kmr), vars = ~ ManEnv),
                      residual = ~ units,
                      weights = "wt.kmr",
                      data = d, maxit = 150, getStartingValues = TRUE)
sv$Init[nrow(sv)] <- 1.0
sv$Fixed[nrow(sv)] <- TRUE
bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ -1 + ManEnv,
                                         random = ~ re(~ exfa(ManEnv, 3):lk(Gen, Ki.kmr), vars = ~ ManEnv),
                                         residual = ~ units,
                                         weights = "wt.kmr",
                                         data = d, maxit = 150,
                                         init = sv$Init, isFixed = sv$Fixed, tolConvLL = tolLL, denseThreshold = 0.1),
                         times = 5)

# MVGK fits ====
rm(list = ls())
suppressMessages(library(cornfruit))
suppressMessages(library(asreml))
suppressMessages(library(kmmr))
suppressMessages(library(Matrix))
suppressMessages(library(microbenchmark))
asreml.options(workspace = "5000mb")

# Loading data:
d <- readRDS("data/DROPS_weights.rds")
D.asr <- readRDS("data/DROPS_ED.rds")[levels(d$Env), levels(d$Env)]
D.kmr <- Matrix(D.asr)
Ki.asr <- solve(readRDS("data/DROPS_K.rds")[levels(d$Gen), levels(d$Gen)])
Ki.kmr <- Matrix(Ki.asr)
attr(Ki.asr, "INVERSE") <- TRUE
d$wt.kmr <- 1 / d$weight

Sys.getenv("MKL_NUM_THREADS")

## MVGK asreml ====
vf <- cornfruit::mvgk(D.asr)
init <- c(rep(var(d$GY) / 20, 28), 0.5, 0.1)
type <- c(rep("V", 28), "R", "V")
cons <- c(rep("P", 28), "U", "P")
bm.asr <- microbenchmark(fit.asr <- asreml(fixed = GY ~ -1 + ManEnv,
                                           random = ~ own(ManEnv, "vf", init, type, cons):vm(Gen, Ki.asr),
                                           residual = ~ units,
                                           weights = weight, family = asr_gaussian(dispersion = 1.0),
                                           data = d, maxit = 150),
                         times = 5)
tolLL <- fit.asr$trace[1, ncol(fit.asr$trace) - 1] - fit.asr$trace[1, ncol(fit.asr$trace) - 2] # 0.00114021503077311

## MVGK kmmr ====
sv <- kmmr(fixed = GY ~ -1 + ManEnv,
           random = ~ re(~ corg(Man):gk(Env, D.kmr):lk(Gen, Ki.kmr), vars = ~ Man:Env),
           residual = ~ units,
           weights = "wt.kmr",
           data = d, maxit = 150,
           getStartingValues = TRUE)
sv$Init[29:31] <- c(0.8, 0.4, 1.0)
sv$Fixed[nrow(sv)] <- TRUE
bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ -1 + ManEnv,
                                         random = ~ re(~ corg(Man):gk(Env, D.kmr):lk(Gen, Ki.kmr), vars = ~ Man:Env),
                                         residual = ~ units,
                                         weights = "wt.kmr",
                                         data = d, maxit = 150,
                                         init = sv$Init, isFixed = sv$Fixed, denseThreshold = 1e-3, tolConvLL = tolLL),
                         times = 5)

## MVGK kmmr GPU ====
sv <- kmmr(fixed = GY ~ -1 + ManEnv,
           random = ~ re(~ corg(Man):gk(Env, D.kmr):lk(Gen, Ki.kmr), vars = ~ Man:Env),
           residual = ~ units,
           weights = "wt.kmr",
           data = d, maxit = 150,
           getStartingValues = TRUE)
sv$Init[29:31] <- c(0.8, 0.4, 1.0)
sv$Fixed[nrow(sv)] <- TRUE
bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ -1 + ManEnv,
                                         random = ~ re(~ corg(Man):gk(Env, D.kmr):lk(Gen, Ki.kmr), vars = ~ Man:Env),
                                         residual = ~ units,
                                         weights = "wt.kmr",
                                         data = d, maxit = 150,
                                         init = sv$Init, isFixed = sv$Fixed, denseThreshold = 1e-3, tolConvLL = tolLL, useGPU = TRUE),
                         times = 5)

# Single stage ====
rm(list = ls())
suppressMessages(library(cornfruit))
suppressMessages(library(asreml))
suppressMessages(library(kmmr))
suppressMessages(library(Matrix))
suppressMessages(library(microbenchmark))
asreml.options(workspace = "5000mb")
K <- readRDS("data/DROPS_K.rds")
d <- readRDS("data/DROPS_singleStage.rds")
d <- droplevels(d[d$type == "CrossValidation",])
d <- droplevels(d[d$Variety_ID %in% rownames(K),])
d <- droplevels(d[d$Experiment != "Deb13R",])
d$Env <- substr(d$Experiment, 1, 5)
d$Man <- substr(d$Experiment, 6, 6)
d$Gen <- d$Variety_ID
d$GY <- d$grain.yield
d$Rep <- d$Replicate
d$Block <- d$block
d$Plot <- d$plotID
d$Row <- d$Row
d$Col <- d$Column
d <- d[, c("Plot", "Row", "Col", "Block", "Rep", "Man", "Env", "Gen", "GY")]
d <- na.omit(d)
d$Rep <- factor(d$Rep, sort(unique(d$Rep)), sort(unique(d$Rep)))
d$Row <- factor(d$Row, sort(unique(d$Row)), sort(unique(d$Row))) # Make sure it's converted to a factor correctly!
d$Col <- factor(d$Col, sort(unique(d$Col)), sort(unique(d$Col))) # Make sure it's converted to a factor correctly!
d$Block <- factor(d$Block, sort(unique(d$Block)), sort(unique(d$Block)))
d$Gen <- as.factor(d$Gen)
d$Env <- as.factor(d$Env)
dss <- droplevels(d[d$Man == "R",])
colnames(dss)[colnames(dss) == "Rep"] <- "compBlock"
colnames(dss)[colnames(dss) == "Block"] <- "incompBlock"
# Then add an actual replicate factor that, for each Env and Gen, uniquely identifies a single record:
dss$Rep <- ave(numeric(nrow(dss)), dss$Env, dss$Gen, FUN = seq_along)
dss$Rep <- factor(dss$Rep, 1:max(dss$Rep))

# Complete blocks are named inconveniently:
dss$compBlock <- as.numeric(dss$compBlock)
dss$compBlock[dss$compBlock == 4] <- 2
dss$compBlock[dss$compBlock == 5] <- 1
dss$compBlock <- as.factor(dss$compBlock)

D.asr <- readRDS("data/DROPS_ED.rds")[levels(dss$Env), levels(dss$Env)]
D.kmr <- Matrix(D.asr)
Ki.asr <- solve(K)
attr(Ki.asr, "INVERSE") <- TRUE
Ki.kmr <- solve(Matrix(K))

Sys.getenv("MKL_NUM_THREADS")

## FA1 asreml ====
m <- 1
bm.asr <- microbenchmark(fit.asr <- asreml(GY ~ -1 + Env,
                                           random = ~ at(Env):idv(compBlock) +
                                             at(Env):idv(compBlock):id(incompBlock) +
                                             at(Env):ar1v(Row):ar1(Col) +
                                             fa(Env, m):vm(Gen, Ki.asr),
                                           residual = ~ dsum(~ units | Env),
                                           data = dss,
                                           maxit = 50),
                         times = 5)
tolLL <- fit.asr$trace[1, ncol(fit.asr$trace) - 1] - fit.asr$trace[1, ncol(fit.asr$trace) - 2] # 0.000235635346143681
u.asr <- fit.asr$coefficients$random[(nrow(fit.asr$coefficients$random) - ((14 + m) * 246) + 1):(nrow(fit.asr$coefficients$random) - (m * 246)), 1]

## FA1 kmmr ====
bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ -1 + Env,
                                         random = ~ atEach("Env", ~ re(~ id(compBlock), vars = ~ 1)) +
                                           atEach("Env", ~ re(~ id(compBlock):id(incompBlock), vars = ~ 1)) +
                                           atEach("Env", ~ re(~ ar1(Row):ar1(Col), vars = ~ 1)) +
                                           re(~ exfa(Env, m):lk(Gen, Ki.kmr), vars = ~ Env),
                                         residual = ~ re(~ id(Env):id(Gen):id(Rep), vars = ~ Env),
                                         data = dss,
                                         stepSize = c(0.1, 0.2, 0.3, 0.4, 0.5),
                                         avInfReg = rep(1.1, 50),
                                         maxit = 50, tolConvLL = tolLL, denseThreshold = 0.1),
                         times = 5)
u.kmr <- fit.kmr$u[(nrow(fit.kmr$u) - ((14 + m) * 246) + 1):(nrow(fit.kmr$u) - (m * 246)), 1]
cor(u.asr, u.kmr)

## FA3 asreml ====
m <- 3
bm.asr <- microbenchmark(fit.asr <- asreml(GY ~ -1 + Env,
                                           random = ~ at(Env):idv(compBlock) +
                                             at(Env):idv(compBlock):id(incompBlock) +
                                             at(Env):ar1v(Row):ar1(Col) +
                                             fa(Env, m):vm(Gen, Ki.asr),
                                           residual = ~ dsum(~ units | Env),
                                           data = dss,
                                           maxit = 50),
                         times = 5)
tolLL <- fit.asr$trace[1, ncol(fit.asr$trace) - 1] - fit.asr$trace[1, ncol(fit.asr$trace) - 2] #  0.00016432644224551
u.asr <- fit.asr$coefficients$random[(nrow(fit.asr$coefficients$random) - ((14 + m) * 246) + 1):(nrow(fit.asr$coefficients$random) - (m * 246)), 1]

## FA3 kmmr ====
bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ -1 + Env,
                                         random = ~ atEach("Env", ~ re(~ id(compBlock), vars = ~ 1)) +
                                           atEach("Env", ~ re(~ id(compBlock):id(incompBlock), vars = ~ 1)) +
                                           atEach("Env", ~ re(~ ar1(Row):ar1(Col), vars = ~ 1)) +
                                           re(~ exfa(Env, m):lk(Gen, Ki.kmr), vars = ~ Env),
                                         residual = ~ re(~ id(Env):id(Gen):id(Rep), vars = ~ Env),
                                         data = dss,
                                         stepSize = c(0.1, 0.2, 0.3, 0.4, 0.5),
                                         avInfReg = rep(1.1, 50),
                                         maxit = 50, tolConvLL = tolLL, denseThreshold = 0.1),
                         times = 5)
u.kmr <- fit.kmr$u[(nrow(fit.kmr$u) - ((14 + m) * 246) + 1):(nrow(fit.kmr$u) - (m * 246)), 1]
cor(u.asr, u.kmr)

# Australian GBLUP ====
rm(list = ls())
suppressMessages(library(cornfruit))
suppressMessages(library(asreml))
suppressMessages(library(kmmr))
suppressMessages(library(sommer))
suppressMessages(library(Matrix))
suppressMessages(library(microbenchmark))
asreml.options(workspace = "10000mb")
d <- read.csv("data/AGT.csv")
d <- d[, c(1, 4)]
colnames(d) <- c("Gen", "GY")
d$Gen <- factor(d$Gen, d$Gen, d$Gen)
Ki.kmr <- readRDS("data/AGT_Ki.rds")
Ki.asr <- as.matrix(Ki.kmr)
Ki.smr <- Ki.kmr
attr(Ki.smr, "inverse") <- TRUE
attr(Ki.asr, "INVERSE") <- TRUE

Sys.getenv("MKL_NUM_THREADS")

## asreml ====
bm.asr <- microbenchmark(fit.asr <- asreml(fixed = GY ~ 1,
                                           random = ~ vm(Gen, Ki.asr),
                                           residual = ~ units,
                                           data = d, maxit = 150),
                         times = 5)
tolLL <- fit.asr$trace[1, ncol(fit.asr$trace) - 1] - fit.asr$trace[1, ncol(fit.asr$trace) - 2] # 0.00141621817601845

## sommer ====
bm.smr <- microbenchmark(fit.smr <- mmes(fixed = GY ~ 1,
                                         random = ~ vsm(ism(Gen), Gu = Ki.smr),
                                         rcov = ~ units,
                                         data = d, nIters = 150, tolParConvLL = tolLL),
                         times = 5)

## kmmr ====
bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ 1,
                                         random = ~ re(~ lk(Gen, Ki.kmr), vars = ~ 1),
                                         residual = ~ units,
                                         data = d, maxit = 150, tolConvLL = tolLL),
                         times = 5)

attr(Ki.kmr, "logDet") <- NULL
bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ 1,
                                         random = ~ re(~ lk(Gen, Ki.kmr), vars = ~ 1),
                                         residual = ~ units,
                                         data = d, maxit = 150, tolConvLL = tolLL),
                         times = 5)

## kmmr GPU ====
bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ 1,
                                         random = ~ re(~ lk(Gen, Ki.kmr), vars = ~ 1),
                                         residual = ~ units,
                                         data = d, maxit = 150, tolConvLL = tolLL, useGPU = TRUE),
                         times = 5)

attr(Ki.kmr, "logDet") <- NULL
bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ 1,
                                         random = ~ re(~ lk(Gen, Ki.kmr), vars = ~ 1),
                                         residual = ~ units,
                                         data = d, maxit = 150, tolConvLL = tolLL, useGPU = TRUE),
                         times = 5)

# DROPS stacking ====
rm(list = ls())
suppressMessages(library(cornfruit))
suppressMessages(library(asreml))
suppressMessages(library(kmmr))
suppressMessages(library(Matrix))
suppressMessages(library(microbenchmark))
asreml.options(workspace = "5000mb")

# Loading data:
d <- readRDS("data/DROPS_weights.rds")
d <- droplevels(d[d$Man == "R",])
D.asr <- readRDS("data/DROPS_ED.rds")[levels(d$Env), levels(d$Env)]
D.kmr <- Matrix(D.asr)
Ki.asr <- solve(readRDS("data/DROPS_K.rds")[levels(d$Gen), levels(d$Gen)])
Ki.kmr <- Matrix(Ki.asr)
attr(Ki.asr, "INVERSE") <- TRUE
d$wt.kmr <- 1 / d$weight

Sys.getenv("MKL_NUM_THREADS")

## MVGK + exfa asreml ====
vf <- cornfruit::mvgk(D.asr)
init <- c(rep(var(d$GY) / 20, 14), 0.1)
type <- c(rep("V", 14), "V")
cons <- c(rep("P", 14), "P")
bm.asr <- microbenchmark(fit.asr <- asreml(fixed = GY ~ -1 + Env,
                                           random = ~ own(Env, "vf", init, type, cons):vm(Gen, Ki.asr) + fa(Env, 2):vm(Gen, Ki.asr),
                                           residual = ~ units,
                                           weights = weight, family = asr_gaussian(dispersion = 1.0),
                                           data = d, maxit = 50), # -630.7962,
                         times = 5)
# tolLL <- fit.asr$trace[1, ncol(fit.asr$trace) - 1] - fit.asr$trace[1, ncol(fit.asr$trace) - 2] # 0.00114021503077311
tolLL <- 0.001

## kmmr stacked CPU ====
sv <- kmmr(fixed = GY ~ -1 + Env,
                random = ~ re(~ stack(~ gk(Env, D.kmr) + fa(Env, 2)):lk(Gen, Ki.kmr), vars = ~ Env),
                residual = ~ units,
                weights = "wt.kmr",
                data = d, maxit = 150,
                getStartingValues = TRUE)
sv$Init[58] <- 1.0
sv$Fixed[nrow(sv)] <- TRUE

bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ -1 + Env,
                                         random = ~ re(~ stack(~ gk(Env, D.kmr) + fa(Env, 2)):lk(Gen, Ki.kmr), vars = ~ Env),
                                         residual = ~ units,
                                         weights = "wt.kmr",
                                         data = d, maxit = 150,
                                         init = sv$Init, isFixed = sv$Fixed, tolConvLL = tolLL,
                                         avInfReg = rep(1.1, 150), denseThreshold = 1e-3, useGPU = FALSE),
                         times = 5)
u.asr <- fit.asr$coefficients$random[1:(2*3276),1]
u.kmr <- fit.kmr$u@x
u.asr.stack <- u.asr[1:3276] + u.asr[3277:(2*3276)]
cor(u.asr.stack, u.kmr)

## kmmr stacked GPU ====
sv <- kmmr(fixed = GY ~ -1 + Env,
           random = ~ re(~ stack(~ gk(Env, D.kmr) + fa(Env, 2)):lk(Gen, Ki.kmr), vars = ~ Env),
           residual = ~ units,
           weights = "wt.kmr",
           data = d, maxit = 150,
           getStartingValues = TRUE)
sv$Init[58] <- 1.0
sv$Fixed[nrow(sv)] <- TRUE

bm.kmr <- microbenchmark(fit.kmr <- kmmr(fixed = GY ~ -1 + Env,
                                         random = ~ re(~ stack(~ gk(Env, D.kmr) + fa(Env, 2)):lk(Gen, Ki.kmr), vars = ~ Env),
                                         residual = ~ units,
                                         weights = "wt.kmr",
                                         data = d, maxit = 150,
                                         init = sv$Init, isFixed = sv$Fixed, tolConvLL = tolLL,
                                         avInfReg = rep(1.1, 150), denseThreshold = 1e-3, useGPU = TRUE),
                         times = 5)
u.asr <- fit.asr$coefficients$random[1:(2*3276),1]
u.kmr <- fit.kmr$u@x
u.asr.stack <- u.asr[1:3276] + u.asr[3277:(2*3276)]
cor(u.asr.stack, u.kmr)

# Single stage DROPS GxExM with Gaussian kernel ====
rm(list = ls())
suppressMessages(library(Matrix))
suppressMessages(library(kmmr))
suppressMessages(library(microbenchmark))

D <- Matrix(readRDS("data/DROPS_ED.rds"))
K <- Matrix(readRDS("data/DROPS_K.rds"))
d <- readRDS("data/DROPS_singleStage.rds")
d <- droplevels(d[d$type == "CrossValidation",])
d <- droplevels(d[d$Variety_ID %in% rownames(K),])
d <- droplevels(d[d$Experiment != "Deb13R",])
d$Env <- substr(d$Experiment, 1, 5)
d$Man <- substr(d$Experiment, 6, 6)
d$Gen <- d$Variety_ID
d$GY <- d$grain.yield
d$Rep <- d$Replicate
d$Block <- d$block
d$Plot <- d$plotID
d$Row <- d$Row
d$Col <- d$Column
d <- d[, c("Plot", "Row", "Col", "Block", "Rep", "Man", "Env", "Gen", "GY")]
d <- na.omit(d)
d$Rep <- factor(d$Rep, sort(unique(d$Rep)), sort(unique(d$Rep)))
d$Row <- factor(d$Row, sort(unique(d$Row)), sort(unique(d$Row))) # Make sure it's converted to a factor correctly!
d$Col <- factor(d$Col, sort(unique(d$Col)), sort(unique(d$Col))) # Make sure it's converted to a factor correctly!
d$Block <- factor(d$Block, sort(unique(d$Block)), sort(unique(d$Block)))
d$Gen <- as.factor(d$Gen)
d$Env <- as.factor(d$Env)
# dss <- droplevels(d[d$Man == "R",])
dss <- d
colnames(dss)[colnames(dss) == "Rep"] <- "compBlock"
colnames(dss)[colnames(dss) == "Block"] <- "incompBlock"
dss$Rep <- ave(numeric(nrow(dss)), dss$Env, dss$Gen, FUN = seq_along)
dss$Rep <- factor(dss$Rep, 1:max(dss$Rep))
dss$compBlock <- as.numeric(dss$compBlock)
dss$compBlock[dss$compBlock == 4] <- 2
dss$compBlock[dss$compBlock == 5] <- 1
dss$compBlock <- as.factor(dss$compBlock)
D <- D[levels(dss$Env), levels(dss$Env)]
K <- K[levels(dss$Gen), levels(dss$Gen)]
Ki <- solve(K)
attr(Ki, "logDet") <- -1 * sum(log(eigen(Ki)$values))
dss$ManEnv <- as.factor(paste(as.character(dss$Man), as.character(dss$Env), sep = ":"))
dss$Man <- as.factor(dss$Man)

Sys.getenv("MKL_NUM_THREADS")

## kmmr CPU ====
bm.kmr <- microbenchmark(fit.kmmr <- kmmr(fixed = GY ~ -1 + ManEnv,
                                          random = ~ atEach("ManEnv", ~ re(~ id(compBlock), vars = ~ 1)) +
                                            atEach("ManEnv", ~ re(~ id(compBlock):id(incompBlock), vars = ~ 1)) +
                                            atEach("ManEnv", ~ re(~ ar1(Row):ar1(Col), vars = ~ 1)) +
                                            re(~ corg(Man):gk(Env, D):lk(Gen, Ki), vars = ~ Man:Env),
                                          residual = ~ re(~ id(Man):id(Env):id(Gen):id(Rep), vars = ~ Man:Env),
                                          data = dss,
                                          stepSize = c(0.1, 0.2, 0.3, 0.4, 0.5),
                                          avInfReg = rep(1.1, 50),
                                          maxit = 50, timing = F, useGPU = F),
                         times = 5)

## kmmr GPU ====
bm.kmr <- microbenchmark(fit.kmmr <- kmmr(fixed = GY ~ -1 + ManEnv,
                                          random = ~ atEach("ManEnv", ~ re(~ id(compBlock), vars = ~ 1)) +
                                            atEach("ManEnv", ~ re(~ id(compBlock):id(incompBlock), vars = ~ 1)) +
                                            atEach("ManEnv", ~ re(~ ar1(Row):ar1(Col), vars = ~ 1)) +
                                            re(~ corg(Man):gk(Env, D):lk(Gen, Ki), vars = ~ Man:Env),
                                          residual = ~ re(~ id(Man):id(Env):id(Gen):id(Rep), vars = ~ Man:Env),
                                          data = dss,
                                          stepSize = c(0.1, 0.2, 0.3, 0.4, 0.5),
                                          avInfReg = rep(1.1, 50),
                                          maxit = 50, timing = F, useGPU = T),
                         times = 5)


