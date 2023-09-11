# **StaGeoMod2022 week 4 -[Generalized Linear Mixed Effect Models Squares (GLMM)]**

## Objectives of this week.

1.	Implement as Generalized Linear Mixed Effect Model that considers the nested structure of the data.
3.  Select the best way to model random effects.
4.	Reduce a Generalized Linear Mixed Effect Model to the minimum adequate number of predictors using a stepwise procedure.

## Generalized Linear Mixed Effect Models (GLMMs) in a nutshell

Where basic statistical methods try to quantify the exact effects of each predictor variable, Ecological and Evolutinary problems often involve random effects, whose purpose is instead to quantify the variation among units. The most familiar types of random effect are the blocks in experiments or observational studies that are replicated across sites or times. Random effects also encompass variation among individuals (when multiple responses are measured per individual, such as survival of multiple offspring or sex ratios of multiple broods), genotypes, species and regions or time periods.

A mixed model, mixed-effects model or mixed error-component model is a statistical model containing both fixed effects and random effects.These models are useful in a wide variety of disciplines in the physical, biological and social sciences. They are particularly useful in settings where repeated measurements are made on the same statistical units (longitudinal study), or where measurements are made on clusters of related statistical units. Because of their advantage in dealing with missing values, mixed effects models are often preferred over more traditional approaches such as repeated measures analysis of variance.

## Why we care?

Up to this point, we have treated all categorical explanatory variables as if they were the same. However, explanatory variables are not all created equal. There are fundamentally different sorts of explanatory variables *fixed effects* and *random effects*.

The short version of the differences between these focuses on each of these explanatory variables ($X_i$) affecting your response variable ($Y$). While **fixed effects** influence only the mean of ($Y$), **random effects** influence only ($Y$) variance. *What does this mean in practical terms?* well, **fixed effects** define the mean response (the mean regression slopes). By comparison, **random effects** govern the variance of the response (you could say the covariance structure of the response variable, and hence their standard errors).

One example of random effects is the hierarchical structure of a sample design that leads to pseudoreplication. To address this, you could concentrate on estimating means of our smaller subset of the data - but that will limit your capacity to draw meaningful conclusions from your data.  Much better to recognize them for what they are, random samples from a much larger population, and to concentrate on their variance - *This is the added variation
caused by differences between the levels of the random effects*.

However, sometimes the added and uneven variability can not be attributed to a random effect. In that case, you would need to model the variability to address the variability problem.

Here is when *Mixed Effects* and * generalised least squares* (GLS) models come to the rescue. The difference is that mixed-effects models allow for nested hierarchical structured datasets and known covariates when modelling the extra variation/heterogeneity. If you use the current data to model the additional variation/heterogeneity, you have a GLS model.
