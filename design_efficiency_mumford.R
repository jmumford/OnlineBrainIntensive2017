## ------------------------------------------------------------------------
# install.packages("neuRosim")
library(neuRosim)

# Let's make a single regressor!
# Run length in seconds
run.length = 60
# TR in seconds
tr = 1
# Onset time in seconds
ons.times = c(10, 20 ,30)
# trial duration in seconds
dur = c(1, 1, 1)
# effect size is for simulating data.  
#  Since we only want a design matrix, 
#  we're setting it to 1 (i.e. don't scale)
eff.size = 1
 
reg = specifydesign(ons.times, dur, run.length, tr,
                    eff.size, conv = "double-gamma")
plot(seq(0, run.length-tr, tr), reg, type = 'l', xlab = "Time (s)",
     ylab = "", main = "Example Regressor")

## ------------------------------------------------------------------------

ntrials.each = 30 # 30 faces and 30 houses

# I will use a while loop to find an order such that 
# there is at most 2 repeats 
keep.looking = 1
while (keep.looking == 1){
  trial.type = sample(rep(c(1,2), each = ntrials.each))
  # The rle function breaks a sequence into lengths and values
  length.repeats = rle(trial.type)$lengths
  # Keep going if the max number of repeats is larger than 2
  keep.looking = max(length.repeats) > 2
}  #end while

dur = rep(1, 2*ntrials.each)
iti = rep(2, 2*ntrials.each)
# Generate onsets using a cumulative sum
# I'm adding a 0 so the first trial starts at 0s
#  Also this will be 1 too long, since the last time
#  is the iti time after the last trial, which doesn't 
#  indicate anything.  This is trimmed off.
ons.all = cumsum(c(0,dur+iti))
ons.all = ons.all[1:(2*ntrials.each)]

run.length = max(ons.all)+9 # 9s beyond last stimulus
# note about run.length.  If it isn't a multiple of TR
# the specifydesign function will fail!!
tr = 1
eff.size = 1
faces = specifydesign(ons.all[trial.type == 1], dur[trial.type == 1],
                      run.length, tr,
                      eff.size, conv = "double-gamma")
houses = specifydesign(ons.all[trial.type == 2], dur[trial.type == 2],
                      run.length, tr,
                      eff.size, conv = "double-gamma")
plot(faces, type = 'l', lwd = 2, col = 'red', xlab = "TRs", 
     ylab = '', ylim = c(min(c(faces, houses)), 1.3))
lines(houses, lwd = 2, col = 'cyan')
legend('topleft', c("faces", "houses"), col = c("red", "cyan"), lwd = c(2,2), bty = 'n')

## ------------------------------------------------------------------------
des.mat = cbind(rep(1, length(faces)), faces, houses)
con = c(0, 1, -1) # faces-houses -> faces > houses

# efficiency
# solve() is how you invert a matrix in R
# R requires %*% for matrix multiplication
1/(t(con)%*%solve(t(des.mat)%*%des.mat)%*%con)

## ------------------------------------------------------------------------

n.loop = 10

# Things that don't change
ntrials.each = 30 
dur = rep(1, 2*ntrials.each)
iti = rep(2, 2*ntrials.each)
ons.all = cumsum(c(0,dur+iti))
ons.all = ons.all[1:(2*ntrials.each)]
run.length = max(ons.all)+9 
tr = 1
eff.size = 1

# Things I'd like to save
eff.val = rep(0, n.loop)
desmats = array(0,c(run.length, 3,  n.loop))

for (i in 1:n.loop){
  keep.looking = 1
  while (keep.looking == 1){
    trial.type = sample(rep(c(1,2), each = ntrials.each))
    length.repeats = rle(trial.type)$lengths 
    keep.looking = max(length.repeats) > 2
  }  #end while

  faces = specifydesign(ons.all[trial.type == 1], dur[trial.type == 1],
                      run.length, tr,
                      eff.size, conv = "double-gamma")
  houses = specifydesign(ons.all[trial.type == 2], dur[trial.type == 2],
                      run.length, tr,
                      eff.size, conv = "double-gamma")
  des.mat = cbind(rep(1, length(faces)), faces, houses)
  con = c(0, 1, -1) # faces-houses -> faces > houses

  eff.val[i] = 1/(t(con)%*%solve(t(des.mat)%*%des.mat)%*%con)
  desmats[,,i] = des.mat
}

# Plot design matrices with best and worst efficiencies
par(mfrow = c(2, 1), mar = c(4, 3, 2, 1))
# best
faces.best = desmats[,2,which(eff.val == max(eff.val))]
houses.best = desmats[,3,which(eff.val == max(eff.val))]
plot(faces.best, type = 'l', lwd = 2, col = 'red', xlab = "TR", 
     ylab = '', ylim = c(min(c(faces.best, houses.best)), 1.3),
     main = "Highest Efficiency")
lines(houses.best, lwd = 2, col = 'cyan')

faces.worst = desmats[,2,which(eff.val == min(eff.val))]
houses.worst = desmats[,3,which(eff.val == min(eff.val))]
plot(faces.worst, type = 'l', lwd = 2, col = 'red', xlab = "TR", 
     ylab = '', ylim = c(min(c(faces.worst, houses.worst)), 1.3),
     main = "Lowest Efficiency")
lines(houses.worst, lwd = 2, col = 'cyan')



## ------------------------------------------------------------------------
# Omit the intercept
diag(solve(cor(des.mat[,2:3])))

## ------------------------------------------------------------------------
# Sample 60 ITI's randomly between 2-4s using a uniform distribution

iti.uni = runif(60, 2, 4)

## ------------------------------------------------------------------------
lambda = 1.5
T = 4
nsamp = 60
shift = 2

m.theory = lambda - T/(exp(1/lambda*T)-1)

R = runif(nsamp)*(1-exp(-T/lambda))
rand.trunc.exp = -log(1-R)*lambda + shift

#theoretical mean
m.theory + shift
#min ITI
min(rand.trunc.exp)
#max ITI
max(rand.trunc.exp)


