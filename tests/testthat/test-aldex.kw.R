library(ALDEx2)

set.seed(1)
data(selex)
group <- c(rep("A", 7), rep("B", 7))
tt <- aldex(selex[1:10,], group, test = "t", mc.samples = 128, bayesEst = FALSE)
gm <- aldex(selex[1:10,], group, test = "kw", mc.samples = 128, bayesEst = FALSe)

test_that("aldex.kw function runs grossly intact", {
  
  expect_equal(
    tt$wi.eBH < .05,
    gm$kw.eBH < .05
  )
  
  expect_equal(
    rownames(tt),
    rownames(gm)
  )
})
