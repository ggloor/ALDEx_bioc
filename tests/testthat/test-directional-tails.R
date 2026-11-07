# Each directional tail includes the observed statistic. The independent
# oracle is stats::wilcox.test, including its paired/tied/normal branches.
test_that("rank tails agree with exhaustive finite null distributions", {
  fast <- if (exists("wilcox.fast")) wilcox.fast else ALDEx2:::wilcox.fast
  for (n1 in 2:4) for (n2 in 2:4) {
    n <- n1+n2
    allocations <- combn(n, n1, simplify=FALSE)
    u <- vapply(allocations, function(a) sum(outer(a, setdiff(seq_len(n), a), ">")), 0)
    for (i in seq_along(allocations)) {
      a <- allocations[[i]]; b <- setdiff(seq_len(n), a)
      data <- rbind(c(a,b), -c(a,b))
      groups <- rep(1:2, c(n1,n2))
      greater <- fast(data, groups, FALSE, "greater")
      less <- fast(data, groups, FALSE, "less")
      expect_equal(greater[1], mean(u >= u[i]), tolerance=1e-12)
      expect_equal(less[1], mean(u <= u[i]), tolerance=1e-12)
      expect_equal(greater, rev(less), tolerance=1e-12)
    }
  }
})

test_that("both fast tails match independent stats tests on every branch", {
  fast <- if (exists("wilcox.fast")) wilcox.fast else ALDEx2:::wilcox.fast
  cases <- list(
    list(a=c(1,1,3), b=c(1,4,5), paired=FALSE),
    list(a=c(1,3,6), b=c(0,1,2), paired=TRUE),
    list(a=c(1,-1,2), b=c(0,0,0), paired=TRUE), # tied absolute differences
    list(a=c(0,1,3), b=c(0,0,0), paired=TRUE), # zero difference
    list(a=seq_len(50), b=51:100, paired=FALSE),
    list(a=seq_len(50), b=rep(0,50), paired=TRUE),
    list(a=c(1,2), b=3:52, paired=FALSE))
  for (case in cases) for (alternative in c("greater", "less")) {
    a <- case$a; b <- case$b
    expected <- suppressWarnings(wilcox.test(a,b,paired=case$paired,
      correct=FALSE,alternative=alternative)$p.value)
    # One feature exercises matrix dimension preservation too.
    actual <- suppressWarnings(fast(matrix(c(a,b),nrow=1),
      rep(1:2,c(length(a),length(b))),case$paired,alternative))
    expect_equal(unname(actual),expected,tolerance=1e-12)
  }
})

test_that("Welch evaluates small tails directly and preserves paired tests", {
  fast <- if (exists("t.fast")) t.fast else ALDEx2:::t.fast
  for (paired in c(FALSE,TRUE)) {
    a <- 1:10; b <- a^2 + 1e6
    data <- rbind(c(a,b), c(b,a))
    actual <- fast(data, rep(1:2,each=10), paired)
    expected <- vapply(seq_len(nrow(data)), function(i)
      t.test(data[i,1:10],data[i,11:20],paired=paired,alternative="less")$p.value,0)
    expect_equal(actual$p.less,expected,tolerance=1e-12)
    expect_gt(actual$p.less[1],0)
    expect_lt(actual$p.less[1],1e-16)
    expect_equal(actual$p,rev(actual$p.less),tolerance=1e-12)
  }
})
