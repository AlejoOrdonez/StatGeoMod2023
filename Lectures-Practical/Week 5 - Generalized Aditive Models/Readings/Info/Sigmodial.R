Chick.1 <- ChickWeight[ChickWeight$Chick == 1, ]
SSlogis(Chick.1$Time, 368, 14, 6)  # response only
Asym <- 368; xmid <- 14; scal <- 6
SSlogis(Chick.1$Time, Asym, xmid, scal) # response and gradient
getInitial(weight ~ SSlogis(Time, Asym, xmid, scal), data = Chick.1)
## Initial values are in fact the converged values
fm1 <- nls(weight ~ SSlogis(Time, Asym, xmid, scal), data = Chick.1)
summary(fm1)

plot(predict(fm1)~ Chick.1$Time)


### Sigmoid function ### create a function to generate sigmoid pattern
sigmoid <- function(x, lower_asymptote, carrying_capacity, growth_rate, time_max) {
    return(lower_asymptote + ((carrying_capacity - lower_asymptote)/(1 + exp(-growth_rate * 
        (x - time_max)))))
}
x <- 1:100
y <- sigmoid(1:100, 1, 50, 0.1, 50) + rnorm(100, 0, 2)
m.s <- nls(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = list(a = min(y), 
    b = max(y), c = 1, d = round(median(x))), trace = TRUE)



plot(y ~ x)
lines(predict(m.s)~x)

a<-((predict(m.s)[-1]-predict(m.s)[-length(x)]))
plot(a[-1]-a[-length(a)])




plot(y ~ x)
lines(predict(m.s)~x)
abline(v=x[which(max(c(a[-1]-a[-length(a)]))==c(a[-1]-a[-length(a)]))-2])

abline(v=x[which(min(c(a[-1]-a[-length(a)]))==c(a[-1]-a[-length(a)]))+2])
