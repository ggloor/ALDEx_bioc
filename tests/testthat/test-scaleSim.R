library(ALDEx2)
set.seed(100)

##Function to create the true abundances via Poisson resampling
create_true_abundances <- function(d, n){
  dd <- length(d)/2
  dat <- as.data.frame(t(sapply(d, function(x) rpois(n, lambda=x))))
  dat <- split(dat, rep(1:2, each=dd))
  dat <- lapply(dat, function(x, dd){`rownames<-`(x, paste0("Taxa", 1:dd))}, dd = dd)
  dat <- lapply(dat, function(x){t(x)})
  dat <- do.call(rbind, dat)
  dat <- as.data.frame(dat)
  dat$Condition <- factor(rep(c("Pre", "Post"), each=n), levels = c("Pre", "Post"))
  dat <- dat[,c(ncol(dat),1:(ncol(dat)-1))]
  rownames(dat) <- NULL
  return(dat)
}

##Function to resample data to an arbitrary size
resample_data <- function(dat, seq.depth){
  ddat <- as.matrix(dat[,-1])/rowSums(as.matrix(dat[,-1]))
  for (i in 1:nrow(dat)){
    dat[i,-1] <- rmultinom(1, size=seq.depth, prob=ddat[i,])
  }
  return(dat)
}


###Setting the data parameters for the simulation
d <- c(4000, 4000, 4000, 4000, 4000, 400,400,400,400,4000,400,500,500,500,400,400,400,400,400,100,400, # Pre
       4000, 4000, 3000, 2000, 4000, 400,400,400,400,4000,400,500,500,500,200,400,400,400,400,100,100) # Post

##Create scale abundances
dat <- create_true_abundances(d, n=50)
##Create resampled data
rdat <- resample_data(dat, seq.depth=5000)

countdata <- t(rdat[,-1,drop=F])
colnames(countdata) <- paste0("n", 1:ncol(countdata))


test_that("aldex2 works without scale samples passed", {
  expect_error(expect_error(aldex(countdata, as.character(rdat$Condition), gamma = NULL, mc.samples = 128))) # expect no error
})

test_that("aldex2 works with scale samples passed", {
  aldex.fit <- aldex(countdata, as.character(rdat$Condition), gamma = .5, mc.samples = 128)
  
  aldex.fit <- aldex.fit[aldex.fit$wi.eBH <= 0.05, ]
  truth <- row.names(aldex.fit)
  expect_true("Taxa21" %in% truth)
})

test_that("aldex2 works with coda scale samples passed", {
  aldex.fit <- aldex(countdata, as.character(rdat$Condition), gamma = 10, mc.samples = 128)
  aldex.fit <- aldex.fit[aldex.fit$wi.eBH <= 0.05, ]
  truth <- row.names(aldex.fit)
  expect_true(length(truth) == 0)
})

## this code was used to set the tolerances
#out.data <- matrix(data=NA, ncol=9, nrow=100)
#
#for(i in 1:100){
#       x <- aldex.clr(selex, conds, verbose=F)
#       x.e <- aldex.effect(x,verbose=F)
#       out.data[i,1:7] <- as.numeric(x.e['P:E:T:E',])
#       out.data[i,8] <- median(x.e$diff.By the way,)
#       out.data[i,9] <- median(x.e$diff.win)
#       print(i)
#}


#########
# unit test
#########

library(ALDEx2)
data(selex)

conds <- c(rep("NS", 7), rep("S", 7))
x <- aldex.clr(selex,conds)
x.e <- aldex.effect(x)

xs <- aldex.clr(selex, conds, gamma=1e-3)
xs.e <- aldex.effect(xs)

xS <- aldex.clr(selex, conds, gamma=1)
xS.e <- aldex.effect(xS)


test_that("scale sim minimally pertubs diff.btw output", {
  
  expect_equal(median(x.e$diff.btw), median(xs.e$diff.btw),
               tolerance=0.1)
  
  expect_equal(median(x.e$diff.win), median(xs.e$diff.win),
               tolerance=0.1)
  
  expect_error(expect_equal(median(x.e$diff.win),
                            median(xS.e$diff.win), tolerance=0.1))
  
})

