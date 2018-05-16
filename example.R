source("amtc.R")

# generate a dataset with 2 changepoints
dt <- 1 : 1e3 * 0.1
dt <- dt + rnorm(length(dt), 0, 5)
dt <- dt + c(rep(0, 300), rep(10, 400), rep(2, 300))
dt <- as.data.frame(dt)
dt$no <- 1 : 1000
dt <- dt[, c(2,1)]
names(dt)[2] <- 'cn'

result <- amtc(dt)
plot(dt, type = 'l')
lines(result@plotDT$fitted, col = 'red')

# generate a dataset with 1 changepoint
dt <- 1 : 1e3 * 0.1
dt <- dt + rnorm(length(dt), 0, 5)
dt <- dt + c(rep(0, 300), rep(10, 700))
dt <- as.data.frame(dt)
dt$no <- 1 : 1000
dt <- dt[, c(2,1)]
names(dt)[2] <- 'cn'
result <- amtc(dt)
plot(dt, type = 'l')
lines(result@plotDT$fitted, col = 'red')

