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
  expect_error(expect_error(aldex(countdata, as.character(rdat$Condition), gamma = NULL, mc.samples = 128, bayesEst = FALSE))) # expect no error
})

test_that("aldex2 works with scale samples passed", {
  aldex.fit <- aldex(countdata, as.character(rdat$Condition), gamma = 1, mc.samples = 128, bayesEst = FALSE)
  
  aldex.fit <- aldex.fit[aldex.fit$wi.eBH <= 0.05, ]
  truth <- row.names(aldex.fit)
  expect_true("Taxa21" %in% truth)
})

test_that("aldex2 works with coda scale samples passed", {
  aldex.fit <- aldex(countdata, as.character(rdat$Condition), gamma = 100, mc.samples = 128, bayesEst = FALSE)
  aldex.fit <- aldex.fit[aldex.fit$wi.eBH <= 0.05, ]
  truth <- row.names(aldex.fit)
  expect_true(length(truth) == 0)
})
