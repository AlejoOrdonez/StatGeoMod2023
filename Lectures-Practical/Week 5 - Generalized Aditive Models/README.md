# **StaGeoMod2022 - week 5 -[Generalized Additive Models (GAMs)]**

## Objectives of this week.

1. Discriminate when? How? and why? Smothers can be used in regression analyses.
2. Implement different types of smothers in `R`.
3. Explore your dataset to avoid common statistical problems when performing a simple/multiple Generalised Additive Model.
4. Implement as simple/multiple Generalised Additive Model in `R`.
5. Determine the fit of a Generalised Additive Model and establish how good a model is.
6. Reduce a model to the minimum adequate number of predictors using a stepwise procedure.

## Generalized Additive Models (GAMs) in a nutshell

In Generalized Linear Models (GLMMS) continuous explanatory variables have been added to models as linear functions, linearized parametric transformations, or through various link functions. In all this cases an explicit or implicit assumption was made about the parametric form of the function to be fitted to the data (whether quadratic, logarithmic, exponential, logistic, reciprocal or whatever). In many cases, however, you have one or more continuous explanatory variables, but you have no a priori reason to choose one particular parametric form over another for describing the shape of the relationship between the response variable and the explanatory variable(s).

Generalized additive models (GAMs) are useful in such cases because they allow you to capture the shape of a relationship between y and x without prejudging the issue by choosing a particular parametric form. Generalized additive models (implemented in `R` by the function `gam()`) extend the range of application of generalized linear models (glm) by allowing non-parametric smoothers in addition to parametric forms, and these can be associated with a range of link functions. All of the error families allowed with glm are available with gam (binomial, poisson, Gamma, etc.).


