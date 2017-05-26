library(ALDEx2)

set.seed(1)
data(selex)
group <- c(rep("A", 7), rep("B", 7))
tt <- aldex(selex[1:10,], group, test = "t", mc.samples = 128)
gm <- aldex(selex[1:10,], group, test = "glm", mc.samples = 128)

test_that("aldex.glm function runs grossly intact", {
  
  expect_equal(
    tt$wi.eBH < .05,
    gm$kw.eBH < .05
  )
  
  expect_equal(
    rownames(tt),
    rownames(gm)
  )
})
