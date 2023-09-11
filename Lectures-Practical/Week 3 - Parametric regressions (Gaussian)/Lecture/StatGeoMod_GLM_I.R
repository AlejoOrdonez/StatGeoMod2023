# Slide 9 - Linear vs non linear modles
dev.new(width=10,height=10)
par(mfrow=c(2,2))
x <- seq(0,10,0.1)y1 <- 2+5*x-0.2*x^2
plot(x,y1,type="l", ylab="y",  main="decelerating\nY~X+X^2",cex.main=2,cex.axis=1.5,cex.lab=2 )


y2 <- 2+5*x-0.4*x^2
plot(x,y2,type="l", ylab="y",  main="Hump\nY~X-X^2",cex.main=2,cex.axis=1.5,cex.lab=2)


y3 <- 2+4*x-0.6*x^2+0.04*x^3
plot(x,y3,type="l", ylab="y",  main="Inflection\nY~X+X^2+X^3",cex.main=2,cex.axis=1.5,cex.lab=2)

y4 <- 2+4*x+2*x^2-0.6*x^3+0.04*x^4
plot(x,y4,type="l", ylab="y", main="local maxima\nY~X+X^2-X^3+X^4",cex.main=2,cex.axis=1.5,cex.lab=2)

dev.new(width=9,height=12)
par(mfrow=c(3,2),oma=rep(0.2,4))
x <- seq(0,10,0.1)
y1 <- exp(x)
plot(x,y1,type="l", ylab="y",  main="Exponential\nY~exp(X)",cex.main=2,cex.axis=1.5,cex.lab=2)

y1 <- log(x)
plot(x,y1,type="l", ylab="y",  main="Logarithmic\nY~log(X)",cex.main=2,cex.axis=1.5,cex.lab=2)

y1 <- 0.5*(1-exp(-2*x))
plot(x,y1,type="l", ylab="y",  main="Asymptotic\nY~a(1-exp(-bx))",cex.main=2,cex.axis=1.5,cex.lab=2)

y1 <- (0.6*X)/(1+0.5*x)
plot(x,y1,type="l", ylab="y",  main="Michaelis–Menten\nY~(ax/(1+bx))",cex.main=2,cex.axis=1.5,cex.lab=2)



x <- seq(0,2*pi,2*pi/100)
y1 <- cos(x)
y2 <- sin(x)
plot(x,y1,type="l", ylab="y",  main="Sine\nY~sine(x)" ,cex.main=2,cex.axis=1.5,cex.lab=2)
plot(x,y2,type="l", ylab="y",  main="Sine\nY~cosine(x)",cex.main=2,cex.axis=1.5,cex.lab=2)





# Slide 12 - R basics - Linear regression
# Build a linear model
Model.1 <- lm(circumference ~ age,
		   data = Orange)


# Slide 13 - R basics - Linear regression Output
# Build a linear model
Model.1 <- lm(circumference ~ age,
		   data = Orange)
Model.1

# Slide 14 - R basics - Linear regression Output
# Build a linear model
Model.1 <- lm(circumference ~ age,
		   data = Orange)
summary(Model.1)

# Slide 15 - R basics - Multiple regression Output
# Build a linear model
Model.2 <- lm(Ozone ~ Solar.R + Wind + Temp,
		   data = airquality)
Model.2

# Slide 16 - R basics - Multiple regression Output
# Build a linear model
Model.2 <- lm(Ozone ~ Solar.R + Wind + Temp,
		   data = airquality)
summary(Model.2)

# Slide 21 - Checking assumptions
# Build a linear model
Model.1 <- lm(circumference ~ age,
		   data = Orange)
par(mfrow=c(2,2))
plot(Model.1)

# Slide 24 - Assumptions - Linearity
plot(circumference ~ age,
	 data = Orange,
	 pch=19, col= "blue",
	 main="circumference Vs age",
	 cex=2,cex.axis=1.5,cex.lab=1.5,cex.sub=1.2,cex.main=2)
Loess.1 <- loess(circumference ~ age,data = Orange)
lines(x=seq(118,1582, by=2),
	  y=predict(Loess.1,
	  			data.frame(age = seq(118,1582, by=2)),
	  			se = F),
	  lwd=2,col="red")

# Slide 25 - Assumptions - Homoscedasticity
# Build a linear model
Model.1 <- lm(circumference ~ age,
		   data = Orange)
plot(Model.1,1,
	pch=19, col= "blue",
	main="Residual Plot",
	cex=2,cex.axis=1.5,cex.lab=1.5,cex.sub=1.2,cex.main=2)
summary(lm(residuals(Model.1)~predict(Model.1)))

# Slide 26 - Assumptions - Homoscedasticity
# Build a linear model
Model.1 <- lm(circumference ~ age,
		      data = Orange)
plot(Model.1,3,
	pch=19, col= "blue",
	main="Stndz. Residual Plot",
	cex=2,cex.axis=1.5,cex.lab=1.5,cex.sub=1.2,cex.main=2)
summary(lm(rstandard(Model.1)~predict(Model.1)))

# Slide 27 - Assumptions - Independence
# Build a linear model
Model.1 <- lm(circumference ~ age,
		   data = Orange)
plot(x= Orange$age,
	 y=rstandard(Model.1),
	 pch=19, col= "blue",
	main="Residual Vs Predictor",
	xlab="Age",ylab="Residual",
	cex=2,cex.axis=1.5,cex.lab=1.5,cex.sub=1.2,cex.main=2)
abline(h=0,lty=2,lwd=1.5)

plot(Model.1,1,
	pch=19, col= "blue",
	main="Residual Plot",
	cex=2,cex.axis=1.5,cex.lab=1.5,cex.sub=1.2,cex.main=2)
	
# Slide 28 - Assumptions - Normality
# Build a linear model
Model.1 <- lm(circumference ~ age,
		      data = Orange)
plot(Model.1,2,
	pch=19, col= "blue",
	main="Quantile - Quantile Plot",
	cex=2,cex.axis=1.5,cex.lab=1.5,cex.sub=1.2,cex.main=2)

# Slide 29 - Assumptions - Leverage
# Build a linear model
Model.1 <- lm(circumference ~ age,
		      data = Orange)
plot(Model.1,5,
     pch=19, col= "blue",
     main="Residual Vs Leverage",
	cex=2,cex.axis=1.5,cex.lab=1.5,cex.sub=1.2,cex.main=2)

# Slide 32 to 35 - Multiple regression assumptions – Multicollinearity
## define colours for each species
my_cols <- c("#00AFBB", "#E7B800", "#FC4E07")  
# Correlation panel
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits=2)
    txt <- paste0("R = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)}
## put histograms on the diagonal
panel.hist <- function(x){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan")}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19, col = my_cols[iris$Species])}
# Create the plots
pairs(iris[,1:4], 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      diag.panel = panel.hist)

# Slide 36 - Multiple regression assumptions – Multicollinearity
Model.3 <- lm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width,
				data=iris)
summary(Model.3)
1-summary(Model.3)$r.square
sapply(1:4,function(i){
			Y <- iris[,i]
			X <- iris[,-i]
			Model.Tol <- lm(Y ~ .,data=X)
			1-summary(Model.Tol)$r.square})

# Slide 42 - Importance of predictors			
Model.2 <- lm(Ozone ~ Solar.R + Wind + Temp,
		   data = airquality)
summary(Model.2)

# Slide 43 - Importance of predictors			
Model.2 <- lm(Ozone ~ Solar.R + Wind + Temp,
		   data = airquality)
Raw.Coefs <- coef(Model.2)[-1]
Sxj <- apply(airquality[,c("Solar.R", "Wind" , "Temp")],2, sd, na.rm=T)
SY <- sd(airquality$Ozone, na.rm=T)
Raw.Coefs*(Sxj/SY)

# Slide 44 - Importance of predictors			
Scle.airquality <- as.data.frame(scale(airquality))
Model.3 <- lm(Ozone ~ Solar.R + Wind + Temp,
		   data = Scle.airquality)
coef(Model.3)
Model.2 <- lm(Ozone ~ Solar.R + Wind + Temp,
		   data = airquality)
Raw.Coefs <- coef(Model.2)[-1]
Sxj <- apply(airquality[,c("Solar.R", "Wind" , "Temp")],2, sd, na.rm=T)
SY <- sd(airquality$Ozone, na.rm=T)
Raw.Coefs*(Sxj/SY)

# Slide 46 - Importance of predictors			
Model.2 <- lm(Ozone ~ Solar.R + Wind + Temp,
		   data = airquality)
summary(Model.2)

# Slide 51 - Do I need this variable?			
Model.2 <- lm(Ozone ~ Solar.R + Wind + Temp,
		   data = airquality)
summary(Model.2)$adj.r.squared
Model.2.Red <- lm(Ozone ~ Solar.R + Wind,
		   			 data = airquality)
summary(Model.2.Red)$adj.r.squared

Model.2 <- lm(Ozone ~ Solar.R + Wind + Temp,
		   data = airquality)
AIC(Model.2)
BIC(Model.2)
Model.2.Red <- lm(Ozone ~ Solar.R + Wind,
		   			 data = airquality)
AIC(Model.2.Red )
BIC(Model.2.Red )


# Slide 52 - Do I need this variable?			
Model.2 <- lm(Ozone ~ Solar.R + Wind + Temp,
		   data = airquality)
Model.2.Red <- lm(Ozone ~ Solar.R + Wind,
		   			 data = airquality)
anova(Model.2.Red,Model.2)

# Slide 52 - Do I need this variable?			
Model.2 <- lm(Ozone ~ Solar.R + Wind + Temp,
		   data = airquality)
Model.2.For <- step(Model.2,
					   direction = "forward",
					   trace = 0)
Model.2.For

Model.2.Bak <- step(Model.2,
					   direction = "backward",
					   trace = 0)
Model.2.Bak

Model.2.Bot <- step(Model.2,
					   direction = "both",
					   trace = 0)
Model.2.Bot