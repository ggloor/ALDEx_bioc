library(ALDEx2)
data(selex)

selex.sub <- selex[1:14, 1:14]
conds <- c(rep("NS", 7), rep("S", 7))
conds3 <- conds
conds3[1:3] <- "D"
mm <- model.matrix(~A, data.frame("A" = conds))
cont <- 1:14

w <- aldex.clr(selex.sub, conds, mc.samples = 2)
x <- aldex.clr(selex.sub, conds3, mc.samples = 2)
y <- aldex.clr(selex.sub, mm, mc.samples = 2)
z <- aldex.clr(selex.sub, cont, mc.samples = 2)

test_that("clr@conds passes along conditions", {
  
  expect_equal(
    w@conds,
    conds
  )
  
  expect_equal(
    x@conds,
    conds3
  )
  
  expect_equal(
    y@conds,
    mm
  )
  
  expect_equal(
    z@conds,
    cont
  )
})
