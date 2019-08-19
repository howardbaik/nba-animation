library(tidyverse)

###########################################################################
###########################################################################
# Load data------
code.dir <- "./code"
data.dir <- "./data"
dat <- read.csv(file=sprintf("%s/2013_11_01_MIA_BKN.csv", data.dir))
###########################################################################
###########################################################################


###########################################################################
###########################################################################
# Plot data for some arbitrary moment------
source(sprintf("%s/constants.R", code.dir)) # loads libraries and constants used throughout code 
source(sprintf("%s/graphics.R", code.dir)) # graphics/plotting functions
# par(mar=c(0, 0, 0, 0)) NOT NEEDED FOR final-result.Rmd
# data.plotter(dat, 1800) NOT NEEDED FOR final-result.Rmd
###########################################################################
###########################################################################

# Transformed Data------
source(sprintf("%s/data_formatting.R", code.dir))
source(sprintf("%s/covariates.R", code.dir))
poss <- possession.indicator(dat) # infer ballcarrier... takes about a minute
tdat <- rearrange.data(dat, poss) # re-shuffle columns by to ballcarrier... (2 min) 
tdat <- offensive.halfcourt(tdat) # transforming to offensive halfcourt
tdat <- offensive.ballcarrier(tdat)
touchID <- get.touchID(tdat)
# `getAllCovars` is inside covariate.R script
covariates <- getAllCovars(tdat) # get covariates... (3 min)
tdat <- data.frame(tdat, touchID=touchID, covariates)


# Since the above section takes a few minutes to complete, 
# its easier to load a pre-computed vrsion of the transformed data set, tdat
load(sprintf("%s/tdat.Rdata", data.dir))

###########################################################################
###########################################################################
# Player similarity adjacency matrix, H------
load(sprintf("%s/playerbases.Rdata", data.dir))
players <- read.csv(sprintf("%s/players2013.csv", data.dir))
head(players)
###########################################################################
###########################################################################

###########################################################################
###########################################################################
# Player's Spatial Occupancy Patterns------
# NOT NEEDED FOR final-result.Rmd

# par(mfrow=c(1,5))
# for(i in 1:5)
#   spatialPlot0(df[i, ], legend=F)
# 
# 
# par(mfrow=c(1,5))
# for(i in 1:5)
#   spatialPlot0(nmf.basis[i, ], legend=F)

###########################################################################
###########################################################################
# Court occupancy distributions------
df.lowrank <- nmf.coef %*% nmf.basis

# NOT NEEDED FOR final-result.Rmd
# par(mfrow=c(1,5))
# for(i in 1:5)
#   spatialPlot0(df.lowrank[i, ], legend=F)


K <- matrix(NA, nrow=nrow(df), ncol=nrow(df)) 

for(i in 1:nrow(K)){
  this.coef <- nmf.coef[i, ] / sum(nmf.coef[i, ])
  K[i, ] <- apply(nmf.coef, 1, function(r) sum((r / sum(r) - this.coef)^2)) 
}

H <- 0 * K
for(i in 1:nrow(H)){
  inds <- order(K[i, ])[1:8 + 1]
  H[i,inds] <- H[inds, i] <- 1 
}



# Check any player’s “neighbors” according to H, we can do (for Al Horford):------
# NOT NEEDED FOR final-result.Rmd

# this.player <- grep("Horford", players$lastname)
# paste(players$firstname, players$lastname)[which(H[this.player, ] == 1)]
# 
# 
# par(mfrow = c(2,5))
# for(i in 1:10)
#   spatialPlot1(take.basis[i, ], legend=F)

###########################################################################
###########################################################################
# Microtransition Model (Plot of acceleration effects as in Figure 4)------
# NOT NEEDED FOR final-result.Rmd

# player.id <- players$player_id[which(players$firstname == "LeBron")] 
# load(sprintf("%s/micros/%s.Rdata", data.dir, player.id))
# component of LeBron James' micro model during ball possession 
# xtable(with.ball$io.x$summary.fixed[, 1:5])

# par(mfrow=c(1,2), mar=c(0,0,0,0))
# vectorPlot(with.ball)
# vectorPlot(without.ball)    

# Code  below estimates the same model parameters for all players on defense
source(sprintf("%s/parameters.R", code.dir)) # loads many modeling functions 
def.micro <- microDefModel(tdat)

# NOT NEEDED FOR final-result.Rmd
# summary(def.micro$mod.x)$coef[, 1:3]


# Macrotransition Entry Models------
load(sprintf("%s/INLA_TAKE.Rdata", data.dir))

this.player <- grep("Bosh", players$lastname)
n.player <- nrow(players)

param.names <- row.names(inla.out$summary.fixed)
n <- nrow(players)
player.params <- matrix(NA, nrow=n, ncol=length(param.names)) 
y.fix <- inla.out$summary.fixed[, "mean"] # fixed effects 
temp <- names(inla.out$summary.random)
basis.inds <- c(which(temp == "p.int"), grep("p.b[0-9][0-9]*", temp)) 
cov.inds <- setdiff(seq(length(inla.out$summary.random)), basis.inds) 

for(pl in 1:n) {
  # add players' random effects to fixed effects
  y.rand <- c(inla.out$summary.random$p.int[pl, "mean"],
              sapply(cov.inds,
                     function(k) inla.out$summary.random[[k]][pl, "mean"]),
              inla.out$summary.random$p.b1[pl + n * (1:n.basis), "mean"])
  player.params[pl, ] <- y.fix + y.rand
}



values <- player.params[this.player, ]
ranks <- apply(player.params, 2, function(col) rank(col)[this.player])


# Chris Bosh's Shot-taking hazard------------
vars <- paste0("b", seq(n.basis))
spat.fixed <- as.numeric(inla.out$summary.fixed["(Intercept)", "mean"] +
                           t(take.basis) %*% inla.out$summary.fixed[vars, "mean"])
spat.random <- as.numeric(inla.out$summary.random$p.int[this.player, "mean"] +
                            t(take.basis) %*% inla.out$summary.random$p.int[this.player + n * (1:n.basis), "mean"])

# NOT NEEDED FOR final-result.Rmd
# par(mfrow=c(1,2), mar=c(1,4,1,6))
# spatialPlot1(spat.fixed + spat.random, axis.args=list(cex.axis=0.75))
# spatialPlot1(spat.random, axis.args=list(cex.axis=0.75))


load(sprintf("%s/INLA_PASS1.Rdata", data.dir))
vars <- paste0("b", seq(n.basis))
spat.fixed <- as.numeric(inla.out$summary.fixed["(Intercept)", "mean"] +
                           t(pass1.basis) %*% inla.out$summary.fixed[vars, "mean"])
spat.random <- as.numeric(inla.out$summary.random$p.int[this.player, "mean"] +
                            t(pass1.basis) %*% inla.out$summary.random$p.int[this.player + n * (1:n.basis), "mean"])

# NOT NEEDED FOR final-result.Rmd
# par(mfrow=c(1,2), mar=c(1,4,1,6))
# spatialPlot2(head(spat.fixed + spat.random, mesh$n),
#              tail(spat.fixed + spat.random, mesh$n),
#              axis.args=list(cex.axis=0.75))
# 


player.id <- players$player_id[grep("Wade", players$lastname)]
load(sprintf("%s/tmats/%s.Rdata", data.dir, player.id))
# names(tmat.ind)



source(sprintf("%s/parameters.R", code.dir))
hyper <- getHyperParams(tdat) # makes sure all parameter inference is loaded 
ev.out <- evLineups(tdat) # coarsened state EVs for each offensive lineup in tdat


lineup.ids <- ev.out$teammates.all[2, ]
this.lineup <- players[match(lineup.ids, players$player_id), ]
# this.lineup[, 2:4]


lineup.states <- paste(rep(this.lineup$lastname, each=14), state_nms) # state names

# EPV curves------
source(sprintf("%s/EPV_calcs.R", code.dir))
source(sprintf("%s/data_formatting.R", code.dir))
source(sprintf("%s/covariates.R", code.dir))

# The below takes a long time
# NOT NEEDED FOR final-result.Rmd
# draw.raw <- multiresDraw(tdat, hyper, def.micro, ev.out, nmic=50, save.positions=F)
# draw <- compressEPV(tdat, draw.raw$fv.epv.list)

load(sprintf("%s/draw.Rdata", data.dir))

# NOT NEEDED FOR final-result.Rmd
# plot(720 - tdat$game_clock[1:100], draw$epv[1:100], xlab="game clock", ylab="EPV")

# every value in draw$epv corresponds to every time value in tdat$game_clock
# * EPV is NA when ball is not in offensive halfcourt with the game clock moving

source(sprintf("%s/EPV_calcs.R", code.dir))
load(sprintf("%s/combined.epv.draws.Rdata", data.dir))

# e.dat IS THE FINAL OUTPUT YOU NEED
e.dat <- combineDatEPV(dat, epv.table)
