test_that("mq does not change the length of the results", {
           d <- mq(1:10, 1, 1:3)    
           expect_that(nrow(d), equals(10))
          })

test_that("mq returns original if d is a function", {
              d <- mq(1:10, 0.4, 1:3, function(d)rep(1,d))
              expect_that( sum(abs(d-(1:10))), equals(0))
          })

test_that("result of mq has the same number of columns as length of q", {
              d <- mq(1:100, c(0.4, 0.6, 0.2), 1:10)
              expect_that(ncol(d), equals(3))
          })

test_that("result of mq has the same number of columns as length of q", {
              d <- mq(1:100, c(0.4, 0.6, 0.2), 1:10)
              expect_that(ncol(d), equals(3))
          })


test_that("mq fails for discount factor length not coinciding with lag structure length", {             
              expect_that(mq(1:10, 1, 2:3, 1), throws_error())
          })

test_that("mq gives warning for non-numeric discount factors", {             
              expect_that(mq(1:10, 1, 2:3, letters[1:2]), gives_warning())
          })

test_that("q8 gives the same results as quantile(type = 8)", {             
              set.seed(123)
              e <- rnorm(5000)
              probs <- seq(0, 1, length.out = 101)
              q1 <- q8(e, probs)
              q2 <- quantile(e, probs, type = 8, names = FALSE)
              expect_that(sum(abs(q1-q2)), equals(0))
          })


test_that("mq gives the same result as mqC without discount factor", {             
              set.seed(1234)
              e <- rnorm(1000)
              probs <- seq(0, 1, length.out = 11)
              m1 <- mqC(e, probs, 1:20, rep(1,20))
              m2 <- mq(e, probs, 1:20)
              expect_that(sum(abs(m1-m2),na.rm=TRUE), equals(0))
          })

test_that("mq gives the same result as mqC with discount factor", {             
              set.seed(1234)
              e <- rnorm(1000)
              d <- runif(20, 0 ,2)
              probs <- seq(0, 1, length.out = 11)
              m1 <- mqC(e, probs, 1:20, d)
              m2 <- mq(e, probs, 1:20, d)
              expect_that(sum(abs(m1-m2),na.rm=TRUE), equals(0))
          })
