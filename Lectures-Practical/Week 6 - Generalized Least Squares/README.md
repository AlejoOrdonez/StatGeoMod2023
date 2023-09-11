# **StaGeoMod2022 week 6 -[Generalized Least Squares (GLSs)]**

## Objectives of this week.

1.	Implement as Generalized Least Squares model that considers:
  a) Unequal variances between observations.
  b) The temporal  correlation between observations. 
  c) Spatial correlation between observations. 
  d) Phylogenetic correlation between observations. 
2.	Evaluate the magnitude of the temporal/spatial/phylogenetic correlation.
3.  Select the best way to model the:
  a) Unequal variances between observations.
  b) temporal/spatial/phylogenetic correlation between observations.
4.	Reduce a Generalized Least Squares model to the minimum adequate number of predictors using a stepwise procedure.

## Generalized Least Squares (GLSs) in a nutshell

The generalized least squares (GLS) estimator of the coefficients of a linear regression is a generalization of the ordinary least squares (OLS) estimator. It is used to deal with situations in which the OLS estimator is not BLUE (best linear unbiased estimator) because one of the main assumptions of the Gauss-Markov theorem, namely that of homoskedasticity and absence of serial correlation, is violated. In such situations, provided that the other assumptions of the Gauss-Markov theorem are satisfied, the GLS estimator is BLUE.

** Ok but why we do this?**

In both ordinary least squares and maximum likelihood approaches to parameter estimation, we made the assumption of constant variance, that is the variance of an observation is the same regardless of the values of the explanatory variables associated with it, and since the explanatory variables determine the mean value of the observation, what we assume is that the variance of the observation is unrelated to the mean.

There are many real situations in which this assumption is inappropriate. In some cases the measurement system used might be a source of variability, and the size of the measurement error is proportional to the measured quantity. Other times this occurs when errors are correlated. Also, when the underlying distribution is continuous, but skewed, such as lognormal, gamma, etc., the variance is not constant, and in many cases variance is a function of the mean.
An important point is that the constant variance is linked to the assumption of normal distribution for the response.

When the assumption of constant variance is not satisfied a possible solution is to transform the data (for example taking log of the response variable and/or the explanatory variables) to achieve constant variance. Another approach is based on generalized or weighted least squares which is an modification of ordinary least squares which takes into account the inequality of variance in the observations. Weighted least squares play an important role in the parameter estimation for generalized linear models.




Text from:

Taboga, Marco (2021). "Generalized least squares", Lectures on probability theory and mathematical statistics. Kindle Direct Publishing. Online appendix. https://www.statlect.com/fundamentals-of-statistics/generalized-least-squares.

http://halweb.uc3m.es/esp/Personal/personas/durban/esp/web/notes/gls.pdf


