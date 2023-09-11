# **StaGeoMod2022 - week 4 -[Generalized Linear Models (GLMs) - ]**

## Objectives of this week.

1.	Implement as simple/multiple Generalized Linear Models in `R`. 
2.	Determine the fit of a Generalized Linear Models and establish how good a model is.
3.	Reduce a Generalized Linear Models to the minimum adequate number of predictors using a stepwise procedure.

## Generalized Linear Models (GLMMS) in a nutshell

So far, most of the analyses described in this course have been based on linear models that assume normally distributed populations of the response variable and the error terms from the fitted models. We deal with a lack of normality in the residuals via transformations on the response variable.

There are situations where transformations are not effective in making errors normal (e.g. when the response variable is categorical) and, in any case, it might be better to model the actual data rather than data that are transformed to meet assumptions. What we need is a technique for modelling that allows other types of distributions
besides normal. Such a technique is called **generalized linear modelling (GLM)**.

Generalized linear models (GLMs) have several characteristics that make them more generally applicable than the general linear models we have considered so far. One of the most important is that least-squares estimation no longer applies, so you must use maximum likelihood methods. As a reminder, Maximum Likelihood Estimation (MLE) methods are based on an iterative search for the parameter's value that makes the observed data most likely to have been observed.

A GLM consists of three components. **First**, is the random component, which is the response variable and its probability distribution. **Second**, is the systematic component, which represents the predictors ($X$ variables) in the model. **Third**, is the link function, which links the random (the response variable or $Y$) and the systematic component (the $X$ variables). The three most common link functions are:

* the identity link ($g_{(\mu)} = \mu$);

* the log link ($g_{(\mu)} = log(\mu)$) used in Poison regressions; 

* the Logit ($g_{(\mu)} = log[\mu/(1 - \mu]$) used in Binomial/logistic regressions.

GLMs are considered parametric models because a probability distribution is specified for the response variable and, therefore, for the error terms from the model.

GLMs are linear models because a linear combination of predictors describes the response variable.

With this fresh in your mind, it is time to start developing a GLMs analysis. In this practical, you will focus on executing a Multiple Logistic regression. **IMPORTANT:** This time around, rather than giving you long explanations, I mainly focus on giving you specific tasks for you to undertake. These tasks indicate the usual analyses process of a dataset with presence-absence data. 

