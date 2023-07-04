library(ALDEx2)

set.seed(1234)
data(selex)
#subset for efficiency
selex <- selex[1201:1240,]
samps <- c(sample(1:7,50,replace = TRUE), sample(8:14, 50, replace = TRUE))
conds <- c(rep("NS", 50), rep("S", 50))
selex <- selex[,samps]
x <- aldex.clr(selex, conds, mc.samples=200, denom="all")
ttest.test <- aldex.ttest(x) 
mm <- model.matrix(~ conds)
x.mm <- x
x.mm@conds <- mm
glm.test <- aldex.glm(x.mm, fdr.method = "BH")

test_that("glm p-values matches t-test", {
  expect_equal(ttest.test$we.ep, glm.test$`condsS:pval`, tol = 1e-2)
  expect_equal(ttest.test$we.eBH, glm.test$`condsS:pval.padj`, tol = 1e-2)
})

test_that("corrected p-values are greater than raw p-values", {
  expect_equal(sum(ttest.test$we.ep > ttest.test$we.eBH), 0)
  expect_equal(sum(ttest.test$wi.ep > ttest.test$wi.eBH), 0)
  expect_equal(sum(glm.test$`condsS:pval` > glm.test$`model.condsS.pval.padj`), 0)
})
