library(ALDEx2)

data(selex)
group <- c(rep("A", 7), rep("B", 7))
clr <- ALDEx2::aldex.clr(selex, group, mc.samples = 128)

# generate the comparison sets from the condition levels
conditions <- group
conditions <- as.factor( conditions )
levels     <- levels( conditions )
levels <- vector( "list", length( levels ) )
names( levels ) <- levels( conditions )
sets <- names(levels)
setAsBinary <- as.numeric(conditions == sets[1])
setA <- which(conditions == sets[1])
setB <- which(conditions == sets[2])
mc.all <- getMonteCarloInstances(clr)
t.input <- sapply(mc.all, function(y){y[, 1]})

test_that("t.fast gives same result as t.test", {
  
  expect_equivalent(
    as.vector(apply(t.input, 1, function(i){
      t.test(x=i[setA],y=i[setB], paired = FALSE)$p.value})),
    ALDEx2:::t.fast(t.input, group, paired = FALSE)
  )
  
  expect_equivalent(
    as.vector(apply(t.input, 1, function(i){
      t.test(x=i[setA],y=i[setB], paired = TRUE)$p.value})),
    ALDEx2:::t.fast(t.input, group, paired = TRUE)
  )
})

test_that("wilcox.fast gives same result as wilcox.test (exact)", {
  
  expect_equivalent(
    as.vector(apply(t.input, 1, function(i){
      wilcox.test(x=i[setA],y=i[setB], paired = FALSE, exact = TRUE)$p.value})),
    ALDEx2:::wilcox.fast(t.input, group, paired = FALSE)
  )
  
  expect_equivalent(
    apply(t.input, 1, function(i){
      wilcox.test(x=i[setA],y=i[setB], paired = TRUE, exact = TRUE)$p.value}),
    ALDEx2:::wilcox.fast(t.input, group, paired = TRUE)
  )
})

data(selex)
dat <- selex
for(i in 1:6){
  fakecol <- data.frame(sample(selex[,1]))
  colnames(fakecol) <- paste0("fake", i)
  dat <- as.data.frame(cbind(dat, fakecol))
}
group <- c(rep("A", 10), rep("B", 10))
clr <- ALDEx2::aldex.clr(dat, group, mc.samples = 128)

# generate the comparison sets from the condition levels
conditions <- group
conditions <- as.factor( conditions )
levels     <- levels( conditions )
levels <- vector( "list", length( levels ) )
names( levels ) <- levels( conditions )
sets <- names(levels)
setAsBinary <- as.numeric(conditions == sets[1])
setA <- which(conditions == sets[1])
setB <- which(conditions == sets[2])
mc.all <- getMonteCarloInstances(clr)
t.input <- sapply(mc.all, function(y){y[, 1]})

test_that("t.fast gives same result as t.test", {
  
  expect_equivalent(
    as.vector(apply(t.input, 1, function(i){
      t.test(x=i[setA],y=i[setB], paired = FALSE)$p.value})),
    ALDEx2:::t.fast(t.input, group, paired = FALSE)
  )
  
  expect_equivalent(
    as.vector(apply(t.input, 1, function(i){
      t.test(x=i[setA],y=i[setB], paired = TRUE)$p.value})),
    ALDEx2:::t.fast(t.input, group, paired = TRUE)
  )
})

test_that("wilcox.fast gives same result as wilcox.test (exact)", {
  
  expect_equivalent(
    as.vector(apply(t.input, 1, function(i){
      wilcox.test(x=i[setA],y=i[setB], paired = FALSE, exact = TRUE)$p.value})),
    ALDEx2:::wilcox.fast(t.input, group, paired = FALSE)
  )
  
  expect_equivalent(
    apply(t.input, 1, function(i){
      wilcox.test(x=i[setA],y=i[setB], paired = TRUE, exact = TRUE)$p.value}),
    ALDEx2:::wilcox.fast(t.input, group, paired = TRUE)
  )
})

data(iris)
dat <- as.data.frame(t(ceiling(iris[1:100, 1:4])))
group <- c(rep("A", 50), rep("B", 50))
clr <- ALDEx2::aldex.clr(dat, group, mc.samples = 128)

# generate the comparison sets from the condition levels
conditions <- group
conditions <- as.factor( conditions )
levels     <- levels( conditions )
levels <- vector( "list", length( levels ) )
names( levels ) <- levels( conditions )
sets <- names(levels)
setAsBinary <- as.numeric(conditions == sets[1])
setA <- which(conditions == sets[1])
setB <- which(conditions == sets[2])
mc.all <- getMonteCarloInstances(clr)
t.input <- sapply(mc.all, function(y){y[, 1]})

test_that("t.fast gives same result as t.test", {
  
  expect_equivalent(
    as.vector(apply(t.input, 1, function(i){
      t.test(x=i[setA],y=i[setB], paired = FALSE)$p.value})),
    ALDEx2:::t.fast(t.input, group, paired = FALSE)
  )
  
  expect_equivalent(
    as.vector(apply(t.input, 1, function(i){
      t.test(x=i[setA],y=i[setB], paired = TRUE)$p.value})),
    ALDEx2:::t.fast(t.input, group, paired = TRUE)
  )
})

test_that("wilcox.fast gives same result as wilcox.test (normal approx.)", {
  
  expect_equivalent(
    as.vector(apply(t.input, 1, function(i){
      wilcox.test(x=i[setA],y=i[setB], paired = FALSE, correct = FALSE)$p.value})),
    ALDEx2:::wilcox.fast(t.input, group, paired = FALSE)
  )
  
  expect_equivalent(
    apply(t.input, 1, function(i){
      wilcox.test(x=i[setA],y=i[setB], paired = TRUE, correct = FALSE)$p.value}),
    ALDEx2:::wilcox.fast(t.input, group, paired = TRUE)
  )
})

t.input[1:2, 1:3] <- t.input[1:2, 1]
t.input[1:2, 51:53] <- t.input[1:2, 51]

test_that("wilcox.fast gives same result as wilcox.test (given ties)", {
  
  expect_equivalent(
    as.vector(apply(t.input, 1, function(i){
      wilcox.test(x=i[setA],y=i[setB], paired = FALSE, correct = FALSE)$p.value})),
    ALDEx2:::wilcox.fast(t.input, group, paired = FALSE)
  )
  
  expect_equivalent(
    apply(t.input, 1, function(i){
      wilcox.test(x=i[setA],y=i[setB], paired = TRUE, correct = FALSE)$p.value}),
    ALDEx2:::wilcox.fast(t.input, group, paired = TRUE)
  )
})
