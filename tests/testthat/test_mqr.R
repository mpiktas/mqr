test_that("mqr gives the same coefficients as lm for models without discount",{
              set.seed(0213)
              res <- numeric(1000)
              e <- rnorm(1000)/5
              res[1:15] <- e[1:15]
              for(i in 16:1000) {
                  q1 <- quantile(res[i-(1:15)], prob=0.5)
                  q2 <- quantile(res[i-(1:15)], prob=1)
                  res[i] <- 1 + 0.2*res[i-1] - 0.3*res[i-12] + 0.4*q1 - 0.2*q2 + e[i]
              }
              yy <- tail(res, 200)
              b <- lm(y ~ mls(y, c(1, 12), 1) + mq(y, 0.5, 1:15) + mq(y, 1, 1:15), data = data.frame(y = yy))
              a <- mqr(y ~ mls(y, c(1, 12), 1) + mq(y, 0.5, 1:15) + mq(y, 1, 1:15), data = data.frame(y = yy))
             
              expect_that(sum(abs(coef(a) - coef(b))), equals(0))
          })

test_that("mqr gives the same summary as lm for models without discount",{
              set.seed(02134)
              res <- numeric(1000)
              e <- rnorm(1000)/5
              res[1:15] <- e[1:15]
              for(i in 16:1000) {
                  q1 <- quantile(res[i-(1:15)], prob=0.5)
                  q2 <- quantile(res[i-(1:15)], prob=1)
                  res[i] <- 1 + 0.2*res[i-1] - 0.3*res[i-12] + 0.4*q1 - 0.2*q2 + e[i]
              }
              yy <- tail(res, 200)
              aa <- summary(mqr(yy ~ mls(yy, c(1, 12), 1) + mq(yy, 0.5, 1:15) + mq(yy, 1, 1:15), data = data.frame(yy = yy)))
              bb <- summary(lm(yy ~ mls(yy, c(1, 12), 1) + mq(yy, 0.5, 1:15) + mq(yy, 1, 1:15), data = data.frame(yy = yy)))
              expect_that(sum(abs(coef(aa) - coef(bb))), equals(0))
          })
