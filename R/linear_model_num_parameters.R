x <- replicate(1, rnorm(100))
y <- 1.2*x[,1]+rnorm(100)
summary(lm.fit <- lm(y~x))
length(lm.fit$coefficients)
length(coef(lm.fit))
logLik(lm.fit)
attributes(logLik(lm.fit))$df
