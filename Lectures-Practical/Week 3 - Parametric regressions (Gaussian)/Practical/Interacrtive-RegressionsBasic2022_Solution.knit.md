---
title: "Multiple Linear Regressions (GLMs)"
documentclass: "report"
output: 
  learnr::tutorial:
  progressive: true
runtime: shiny_prerendered
---



## **The Set-up**

*Learning objectives.*
1. Explore your dataset to avoid common statistical problems when performing a simple/multiple linear regression.
2.	Implement as simple/multiple linear regression in `R`.  
3.	Discriminate between raw and standardised regressions coefficients. 
4.	Determine the fit of a linear regression and establish how good a model is.  
5.	Reduce a model to the minimum adequate number of predictors using a stepwise procedure.   

## **Environmental factors driving Denmark's biodiversity patterns.**

**About the data you will use.**

The data to be used in this week practical is a combination of the biodiversity data  from the [Biowide](https://ecos.au.dk/forskningraadgivning/temasider/biowide) project, which aimed at collecting biodiversity data from 130 terrestrial sampling sites across Denmark. The dataset contains data for Plants (vascular plants and bryophytes), fungi (lichens and macrofungi), gastropods, and arthropods. You can get the [biodiversity data](https://www.gbif.org/dataset/cd4eb3ec-0a18-42b4-bda4-155716ddd7b1#description) from the Danish Biodiversity Information Facility. For further information on the project check the publication.

Although environmental variables were collected as part of Biowide, this data is not available for our use. The data on environmental variables included in the table you will use comes from multiple Geospatial datasets I complied for this course.

Here you will focus on the flowing variables:

**1) decimalLongitude**: East-West position in decimal-degrees of a sample site coming from [Biowide](https://ecos.au.dk/forskningraadgivning/temasider/biowide).

**2) decimalLatitude**: North-South position in decimal-degrees of a sample site coming from [Biowide](https://ecos.au.dk/forskningraadgivning/temasider/biowide). 

**3) Log.S.AllGrp**: Logarithm of the total species richness of a sample site (sum of the unique species of all sampled groups in a site) coming from [Biowide](https://ecos.au.dk/forskningraadgivning/temasider/biowide).

**4) landsdele**: Regional political unit of Denmark coming form the Database of Global Administrative Areas [GADAM](https://gadm.org).

**5) Geology**: Geological formation coming from the [Denamrk's geological service](https://eng.geus.dk/products-services-facilities/data-and-maps/maps-of-denmark).

**6) NatDensBasMp**: Density of natural land cover types as defined by [BASEMAP](https://envs.au.dk/en/research-areas/society-environment-and-resources/land-use-and-gis/basemap).

**7) Dis2Cost**: Distance in km of a sampled site to the coast estimated by me using Denmark high resolution contour map coming form the Database of Global Administrative Areas [GADAM](https://gadm.org).

**8) HII**: [Human Influence Index](https://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-influence-index-geographic), a measurement of the anthropogenic impacts on the environment.

**9) Slope30mAgg**: the maximum rate of change in value from that cell to its neighbours estimated based on the [SRTM-30m DEM](https://www2.jpl.nasa.gov/srtm/).

**10) PrecSea**: Seasonality of precipitation measured as the variability of monthly values (using the CV) coming the [CHELSA-climatologies](https://chelsa-climate.org) .


## 1. Load the Data

<div class="alert alert-info">
**Your task:**
Load the data in *Biowide_AllSppRich.csv*, located in the files of the project

Use `read.csv` to load the file - The file is comma-separated

Save the object as an object named `DK.Biodiv`.
</div >

<div class="tutorial-exercise" data-label="LoadData" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Read the file and save it as an object named `DK.Biodiv`
DK.Biodiv <-  -----------(-----------)
# print the Five first records
---------(------)
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="LoadData-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Read the file and save it as an object named `DK.Biodiv`
DK.Biodiv <-  read.csv("Data/Biowide_AllSppRich.csv")
# print the Five first records
head(DK.Biodiv)
```

</div>


## 2. Data Exploration: oultiers in X and Y

<div class="alert alert-info">
**Your task:**

Using the object `DK.Biodiv`, you will now build a box plot for each continuous predictor variables adding a "title" to each plot using the `main` argument.

To be clear the response variable is  *Richness of all sampled species* [`Log.S.AllGrp`].

The predictors are:
1) *Region* [`landsdele`]
2) *Parental material* [`Geology`]
3) *Nature Density* [`NatDensBasMp`]
4) *Distance to coast* [`Dis2Cost`]
5) *Human Impact Index* [`HII`]
6) *Slope* [`Slope30mAgg`]
7) *Precipitation Seasonality* [`PrecSea`]
</div>

<div class="tutorial-exercise" data-label="Outlier1" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="15">

```text
# Set a plotting space with 6 regions (two columns and three rows). For this, use the function par()
par(mfrow = c(3,2))

# Plot a box-plot for Nature Density
boxplot(-----------,
        main = "Nature Density")
# Plot a box-plot for Distance to coast
-----------(-----------,
        ----------- = "-----------")
# Plot a box-plot for Human Impact Index
-----------(-----------,
        ----------- = "-----------")
# Plot a box-plot for Slope
-----------(-----------,
        ----------- = "-----------")
# Plot a box-plot for Precipitation Seasonality
-----------(-----------,
        ----------- = "-----------")
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="Outlier1-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Set a plotting space with 6 regions (two columns and three rows). For this, use the function par()
par(mfrow = c(3,2),mar=rep(2,4))

# Plot a box-plot for Nature Density
boxplot(DK.Biodiv$NatDensBasMp,
        main = "Nature Density")
# Plot a box-plot for Distance to coast
boxplot(DK.Biodiv$Dis2Cost,
        main = "Distance to coast")
# Plot a box-plot for Human Impact Index
boxplot(DK.Biodiv$HII,
        main = "Human Impact Index")
# Plot a box-plot for Slope
boxplot(DK.Biodiv$Slope30mAgg,
        main = "Slope")
# Plot a box-plot for Precipitation Seasonality
boxplot(DK.Biodiv$PrecSea,
        main = "Precipitation Seasonality")
```

</div>


```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="OutlierQ1" class="tutorial-question panel-body">
<div id="OutlierQ1-answer_container" class="shiny-html-output"></div>
<div id="OutlierQ1-message_container" class="shiny-html-output"></div>
<div id="OutlierQ1-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```

## 3. Data Exploration: Is there homogeneity of variances - 1?


<div class="alert alert-info">
**Your task:**

Here you will evaluate if the variance in the richness of all sampled species in the **Biowide** project are homogeneous across all predictors. **Remember that the response variable here is `Log.S.AllGrp`**.

For this, you will assess this assumption for each of the five predictors plus Latitude and longitude **individually**. 

You should use a `for()` loop to determine homoscedasticity of Log-species richness  across predictors.

</div>

<div class="tutorial-exercise" data-label="HomVar2" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="25">

```text
# Build 4x3 plotting space filled row-wise 
par(mfrow = c(4,2),
    mar = rep(4,4))

# Create a Vector named PredNames with the predictors names
PredNames <- ----------

# Loop using i as an iterator and cycle through PredNames
for (i in ----------){ 
# Create a new subset data.frame named DK.Div.Tmp that only contains ONLY the log-richness and the evaluated predictors as variables
  DK.Div.Tmp <- ----------[,c("----------",--)]
  
# Build a regression model, using the lm() function, and store it as an object named Lm.Mod 
  Lm.Mod <- lm(Log.S.AllGrp ~ ., # this is a special format for formulas
               data = ----------) 
  
# Save the residuals into an object named Lm.Resid
  Lm.Resid <- ----------(----------)

# Save the fitted/predicted values into an object named Lm.Fitt
  Lm.Fitt <- ----------(----------)

  # Plot the Residuals vs Predicted values
  plot(---------- ~ ----------,
       pch=19, # This argument as defined sets the plotting point to a filled dot.
       xlab = "Fitted values", # This argument define the text to be added as the x-axis label
       ylab = "Residuals", # This argument define the text to be added as the y-axis label
       main = i) # This argument define the text to be added figure title
# Add a horizontal line using the abline() function at zero for reference
  ----------( h = 0)
}
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="HomVar2-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Build 4x3 plotting space filled row-wise 
par(mfrow = c(4,2),
    mar = rep(2,4))

# Create a Vector named PredNames with the predictors names
PredNames <- names(DK.Biodiv)[-c(1,4:6)]

# Loop using i as an iterator and cycle through PredNames
for (i in PredNames){ 
# Create a new subset data.frame named DK.Div.Tmp that only contains ONLY the log-richness and the evaluated predictors as variables
  DK.Div.Tmp <- DK.Biodiv[,c("Log.S.AllGrp",i)]
  
# Build a regression model, using the lm() function, and store it as an object named Lm.Mod 
  Lm.Mod <- lm(Log.S.AllGrp ~ ., 
               data = DK.Div.Tmp) 
  
# Save the residuals into an object named Lm.Resid
  Lm.Resid <- residuals(Lm.Mod)

# Save the fitted/predicted values into an object named Lm.Fitt
  Lm.Fitt <- predict(Lm.Mod)

  # Plot the Residuals vs Predicted values
  plot(Lm.Resid ~ Lm.Fitt,
       pch=19, # This argument as defined sets the plotting point to a filled dot.
       xlab = "Fitted values", # This argument define the text to be added as the x-axis label
       ylab = "Residuals", # This argument define the text to be added as the y-axis label
       main = i) # This argument define the text to be added figure title
# Add a horizontal line using the abline() function at zero for reference
  abline(h=0)
}
```

</div>

```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="HomVarQ1" class="tutorial-question panel-body">
<div id="HomVarQ1-answer_container" class="shiny-html-output"></div>
<div id="HomVarQ1-message_container" class="shiny-html-output"></div>
<div id="HomVarQ1-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```

## 4. Data Exploration: Is there homogeneity of variances - 2?

<div class="alert alert-info">
**Your task:**

Now you will evaluate the variance in the log-richness for all sampled species is homogeneous across all predictors **simultaneously**.

For this, you will build a multiple regression where all seven(7) continuous predictors are combined into a single additive model **with no interactions**.

**Remember that your response variable is `Log.S.AllGrp`**.
</div>

<div class="tutorial-exercise" data-label="HomVar1" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="20">

```text
# Build a regression model using the lm() function where all six variables predict the log-richness for all sampled species. # save the model as an object named FullMod
FullMod <- lm(---------- ~ ---------- + ---------- + ---------- + ---------- + ---------- + ---------- + ----------,
              data = ----------)
# Call the model
FullMod

# Now extract the residuals and save them into an object named Lm.Resid. For this, use the residuals() function
Lm.Resid <- ----------(----------)

# Now extract the fitted values and save them into an object named Lm.Fitt.  For this, use the predict() function
Lm.Fitt <- ----------(----------)

# Plot the Residuals vs Predicted values
  plot(---------- ~ ----------,
       pch=19, # This argument, as defined, sets the plotting point to a filled dot.
       xlab = "Fitted values", # This argument defines the text to be added as the x-axis label
       ylab = "Residuals", # This argument define the text to be added as the y-axis label
       main = "All variables", # This argument define the text to be added figure title
      ylim = range(Lm.Resid)) 
# Add a horizontal line using the abline() function at zero for reference  
  ----------(---------- = ----------)
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="HomVar1-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Build a regression model using the lm() function where all six variables predict the log-richness for all sampled species. # save the model as an object named FullMod
FullMod <- lm(Log.S.AllGrp ~ NatDensBasMp + Dis2Cost + HII + Slope30mAgg + PrecSea + decimalLongitude + decimalLatitude,
              data = DK.Biodiv)
# Call the model
FullMod

# Now extract the residuals and save them into an object named Lm.Resid. For this, use the residuals() function
Lm.Resid <- residuals(FullMod)

# Now extract the fitted values and save them into an object named Lm.Fitt.  For this, use the predict() function
Lm.Fitt <- predict(FullMod)

# Plot the Residuals vs Predicted values
  plot(Lm.Resid ~ Lm.Fitt,
       pch=19, # This argument, as defined, sets the plotting point to a filled dot.
       xlab = "Fitted values", # This argument defines the text to be added as the x-axis label
       ylab = "Residuals", # This argument define the text to be added as the y-axis label
       main = "All variables", # This argument define the text to be added figure title
      ylim = range(Lm.Resid)) 
# Add a horizontal line using the abline() function at zero for reference  
  abline(h = 0)
```

</div>


```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="HomVarQ2" class="tutorial-question panel-body">
<div id="HomVarQ2-answer_container" class="shiny-html-output"></div>
<div id="HomVarQ2-message_container" class="shiny-html-output"></div>
<div id="HomVarQ2-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```

## 5a. Data Exploration: Are the residuals normally distributed?

<div class="alert alert-info">
**Your task:**

Using a simple multivariate additive model (that is and additive combination of the used predictors) build in the section above (that is the `Lm.Resid` object), you will now assess if the residuals are normally distributed. For this, you will:

* Generate a histogram (using the `hist()` function) of the **frequency** of the residuals from the simple multivariate additive model.
</div>
<div class="tutorial-exercise" data-label="NormRes2" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Use the hist() function to display the distribution of the residuals. Here do not plot the count but the frequency.
hist(----------,
     freq = ----)

# Estimate the density of the residuals using the function `density()`
Den.Lm.Resid <- ----------(----------)

# Plot the density estimates of the residuals over the histogram using the function `lines()`
----------(----------)
#
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="NormRes2-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Use the hist() function to display the distribution of the residuals. Here do not plot the count but the frequency.
hist(Lm.Resid,
     freq = F)

# Estimate the density of the residuals using the function `density()`
Den.Lm.Resid <- density(Lm.Resid)

# Plot the density plot over the histogram

lines(Den.Lm.Resid)
#
```

</div>


```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="NormResQ1" class="tutorial-question panel-body">
<div id="NormResQ1-answer_container" class="shiny-html-output"></div>
<div id="NormResQ1-message_container" class="shiny-html-output"></div>
<div id="NormResQ1-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```

## 5b. Data Exploration: Are the residuals normally distributed?

<div class="alert alert-info">
**Your task:**
Using a simple multivariate additive model (that is and additive combination of the used predictors) build in the section above (that is the `Lm.Resid` object), you will now assess if the residuals are normally distributed. For this, you will:

* Plot a q-q plot of the residuals, using the `qqnorm()` and `qqline()` functions.

</div>

<div class="tutorial-exercise" data-label="NormRes3" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# use the qqnorm() on the object containing the residuals to make a q-q plot
----------(----------)
# use the qqline in the object containing the residuals to add the expectation line
----------(----------)
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="NormRes3-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# use the qqnorm() and qqline() functions to plot a qqplot
qqnorm(Lm.Resid) # use the qqnorm() on the object containing the residuals to make a q-q plot
qqline(Lm.Resid) # use the qqline in the object containing the residuals to add the expectation line
```

</div>

```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="NormResQ2" class="tutorial-question panel-body">
<div id="NormResQ2-answer_container" class="shiny-html-output"></div>
<div id="NormResQ2-message_container" class="shiny-html-output"></div>
<div id="NormResQ2-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```

## 5c. Data Exploration: Are the residuals normally distributed?
<div class="alert alert-info">
**Your task:**

Using a simple multivariate additive model (that is and additive combination of the used predictors) build in the section above (that is the `Lm.Resid` object), you will now assess if the residuals are normally distributed. For this, you will:

* Use the Shapiro–Wilk test (executed using the `shapiro.test()` function) to assess the normality of the residuals.
</div>

<div class="tutorial-exercise" data-label="NormRes4" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Implement the shapiro.test() function to assess normality
----------(----------)
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="NormRes4-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Implement the shapiro.test() function to assess normality of the residuals. Remember that a non-significant test (p>0.05) means that the observed and expected distribution (in this case, normal) are the same.
shapiro.test(Lm.Resid)
```

</div>

```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="NormResQ3" class="tutorial-question panel-body">
<div id="NormResQ3-answer_container" class="shiny-html-output"></div>
<div id="NormResQ3-message_container" class="shiny-html-output"></div>
<div id="NormResQ3-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```

## 6. Data Exploration: Is there collinearity among the covariates?

<div class="alert alert-info">
**Your task:**

Estimate the tolerance and VIF for all six predictors in the `DK.Biodiv` object. Use a Loop to estimate the tolerances for all variables.

</div>

<div class="tutorial-exercise" data-label="Collin3" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="20">

```text
# Create a vector with only NA values named Tol.Summ to store the tolerance values
Tol.Summ <- rep(x = NA,
                times = length(PredNames))
# give names to each place in the vector using PredNames
----------(----------) <- PredNames

# Use a for loop across predictors (remember the names are in PredNames) to estimate the tolerance.
for (i in PredNames){
# Create new data.frame named Tol.DF was the first variable tested.
  Tol.DF <- data.frame(Pred = DK.Biodiv[,---], # the data of the predictor variable being tested
                       DK.Biodiv[,----------])# Add the other predictor variables here.   # Build an lm model to predict the predictor of interest as a function of all other predictors. Save this model as Tol.LM.
    Tol.LM <- lm(Pred~.,
                 data = Tol.DF) 
# To estimate the tolerance for the predictor of interest, the first step is to extract the R2 from the regression model and store it as an object name Tol.R2
    Tol.R2 <- ----------(----------)$---------- #

# estimate the tolerance (1-R2) and save it in the corresponding position of the Tol.Summ you created before the loop. 
    Tol.Summ[i] <- 1 - Tol.R2
}    

# now Call the vector Tol.Summ you created before the loop 
----------
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="Collin3-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Create a vector with only NA values named Tol.Summ to store the tolerance values
Tol.Summ <- rep(x = NA,
                times = length(PredNames))
# give names to each place in the vector using PredNames
names(Tol.Summ) <- PredNames

# Use a for loop across predictors (remember the names are in PredNames) to estimate the tolerance.
for (i in PredNames){
# Create new data.frame named Tol.DF was the first variable tested.
  Tol.DF <- data.frame(Pred = DK.Biodiv[,i], # the data of the predictor variable being tested should be included here  - call it using the i iterator.
                       DK.Biodiv[,PredNames[!PredNames%in%i]])# Add the other predictor variables here. For this, use the logical test !PredNames%in%i (each select all predictors variable EXECEPT i). 
  
  # Build an lm model to predict the predictor of interest as a function of all other predictors. Save this model as Tol.LM.
    Tol.LM <- lm(Pred~.,
                 data = Tol.DF) 
# To estimate the tolerance for the predictor of interest, the first step is to extract the R2 from the regression model and store it as an object name Tol.R2
    Tol.R2 <- summary(Tol.LM)$r.squared #

# estimate the tolerance (1-R2) and save it in the corresponding position of the Tol.Summ you created before the loop. 
    Tol.Summ[i] <- 1 - Tol.R2
}    

# now Call the vector Tol.Summ you created before the loop 
Tol.Summ
```

</div>


```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="CollinQ1" class="tutorial-question panel-body">
<div id="CollinQ1-answer_container" class="shiny-html-output"></div>
<div id="CollinQ1-message_container" class="shiny-html-output"></div>
<div id="CollinQ1-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```

## 7. Data Exploration: What are the relationships between `Y` and `X` variables -1?

<div class="alert alert-info">
**Your task:**

Make multi-panel scatter-plots showing the relation between log-richness for all sampled species and each predictors.

For each panel add the regression line (using the `abline()` function) and the Pearson correlation coefficient (using the `legend()` function).

</div>

<div class="tutorial-exercise" data-label="BivarRel" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="20">

```text
# Build 3x2 plotting space filled row-wise 
par(mfrow = c(4,2), 
    mar=c(4,4,2,2))

# Now lets loop across predictors to estimate the Tolerance.
for (i in PredNames){

  # Subset DK.Biodiv to create a new data.frame named Temp.DF that has Log.S.AllGrp and the predictor evaluated
  Temp.DF <- ----------[,c(----------,"----------")]

  # Plot the bivariate relation
  plot(x = ----------, # Call the predictor
       y = ----------, # Call the response (Log.S.AllGrp)
       pch = 19, # Set the points to a filled dot
       xlab = i, # Set the Name for the X-axis
       ylab = "Log.S.AllGrp") # Set the Name for the Y-axis

# Add the regression line - using the abline function
# First create a regression object named Reg.Temp that has the regression between Log.S.AllGrp and the predictor
Reg.Temp <- lm(---------- ~ ----------)

# Using the abline function plot the object with the regression object
abline(----------, # The regression object
       col="red", # Set the line colour to "red"
       lwd=1.5) # make the plotted line thicker

# Add the the Pearson correlation - make the text boldface if the relation is significant
legend("topleft",
       legend = ----------,
       text.font = ----------)
}
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="BivarRel-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Build 3x2 plotting space filled row-wise 
par(mfrow = c(4,2), 
    mar=c(4,4,2,2))

# Now lets loop across predictors to estimate the Tolerance.
for (i in PredNames){

  # Subset DK.Biodiv to create a new data.frame named Temp.DF that has Log.S.AllGrp and the predictor evaluated
  Temp.DF <- DK.Biodiv[,c(i,"Log.S.AllGrp")]

  # Plot the bivariate relation
  plot(x = Temp.DF[,i], # Call the predictor
       y = Temp.DF$Log.S.AllGrp, # Call the response (Log.S.AllGrp)
       pch = 19, # Set the points to a filled dot
       xlab = i, # Set the Name for the X-axis
       ylab = "Log.S.AllGrp") # Set the Name for the Y-axis


# Add the regression line - using the abline function
# First create a regression object named Reg.Temp that has the regression between Log.S.AllGrp and the predictor
Reg.Temp <- lm(DK.Biodiv[,"Log.S.AllGrp"] ~ DK.Biodiv[,i])

# Using the abline function plot the object with the regression object
abline(Reg.Temp, # The regression object
       col="red", # Set the line colour to "red"
       lwd=1.5) # make the plotted line thicker

# Add the the Pearson correlation - make the text boldface if the relation is significant
legend("topleft",
       legend = paste0("r = ",
                       round(cor(na.omit(Temp.DF))[2,1],3)
                       ),
       text.font = ifelse(cor.test(na.omit(Temp.DF)[,1],na.omit(Temp.DF)[,2])$p.value<0.05,2,1))
}
```

</div>


```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="BivarRelQ1" class="tutorial-question panel-body">
<div id="BivarRelQ1-answer_container" class="shiny-html-output"></div>
<div id="BivarRelQ1-message_container" class="shiny-html-output"></div>
<div id="BivarRelQ1-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```

## 8. Data Exploration: Are observations of the response variable independent?

<div class="alert alert-info">
**Your task:**
Now you will evaluate if the independence of the Observations. For this, plot the residuals (`Lm.Resid`) vs. the response variable (`Log.S.AllGrp`).
</div>

<div class="tutorial-exercise" data-label="Independ1" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Plot the Residuals vs response variable
  plot(---------- ~ ----------,
       data = ----------,
       pch=19, # This argument, as defined, sets the plotting point to a filled dot.
       xlab = "Log.S.AllGrp", # This argument defines the text to be added as the x-axis label
       ylab = "Residuals", # This argument define the text to be added as the y-axis label
       main = "All variables", # This argument define the text to be added figure title
      ylim = ----------(----------)) 
# Add a horizontal line using the abline() function at zero for reference
  ----------(----------)
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="Independ1-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Plot the Residuals vs response variable
  plot(Lm.Resid ~ Log.S.AllGrp,
       data = DK.Biodiv[as.numeric(names(Lm.Resid)),],
       pch=19, # This argument, as defined, sets the plotting point to a filled dot.
       xlab = "Log.S.AllGrp", # This argument defines the text to be added as the x-axis label
       ylab = "Residuals", # This argument define the text to be added as the y-axis label
       main = "All variables", # This argument define the text to be added figure title
      ylim = c(-max(abs(Lm.Resid)), # a new argument you will use here is ylim to set the lower and upper limits of the y-axes. For this, you will specify a vector of two-element with the lower and upper limits. This will allow you to make the y-axis symmetric.
               max(abs(Lm.Resid)))) 
# Add a horizontal line using the abline() function at zero for reference
  abline(h=0)
```

</div>


```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="IndependQ1" class="tutorial-question panel-body">
<div id="IndependQ1-answer_container" class="shiny-html-output"></div>
<div id="IndependQ1-message_container" class="shiny-html-output"></div>
<div id="IndependQ1-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```


## 9. Multiple Linear Regressions - first implementation

<div class="alert alert-info">
**Your task:**
You have a regression object (`FullMod`). Plot the regression object to see what is shown.
</div>

<div class="tutorial-exercise" data-label="Reg2" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Build 2x2 plotting space filled row-wise 
par(mfrow = c(2,2), 
    mar=c(4,4,2,2))
# use the plot function on the lm object named FullMod
----------(----------)
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="Reg2-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Build 2x2 plotting space filled row-wise 
par(mfrow = c(2,2), # you need to specify two values here the number of rows and the number of columns
    mar=c(4,4,2,2)) # here you define the 'margins" - blank space between plotting areas)
# use the plot function on the lm object named FullMod
plot(FullMod)
```

</div>

## 10. Multiple Linear Regressions - standardised regression coefficients

<div class="alert alert-info">
**Your task:**
Create a new `data.frame` including `Log.S.AllGrp` and all six predictors, and use the `scale()` function to standardise both the response and the predictor variables.
</div>

<div class="tutorial-exercise" data-label="StdzPred1" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# For ease, generate an object named Stdz.DK.Biodiv where you will store `Log.S.AllGrp` and the used of predictors
Stdz.DK.Biodiv <- ----------[,c("----------", # the name of the response
                            ----------)]# the name of the predictors
# Use the scale() function on Stdz.DK.Biodiv to scale the predictors
Stdz.DK.Biodiv <- ----------(----------)

# Print the Five first rows of Stdz.DK.Biodiv
----------(----------)
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>

<div class="tutorial-exercise-support" data-label="StdzPred1-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# For ease, generate an object named Stdz.DK.Biodiv where you will store `Log.S.AllGrp` and the used of predictors
Stdz.DK.Biodiv <- DK.Biodiv[,c("Log.S.AllGrp", # the name of the response
                            PredNames)]# the name of the predictors
# Use the scale() function on Stdz.DK.Biodiv to scale the predictors
Stdz.DK.Biodiv <- scale(Stdz.DK.Biodiv)

# Print the Five first rows of Stdz.DK.Biodiv
head(Stdz.DK.Biodiv)
```

</div>
<div class="alert alert-info">
**Your task:**

Estimate the standardised regression coefficients using the standardised `data.frame` you just created (`Stdz.DK.Biodiv`).

</div>

<div class="tutorial-exercise" data-label="StdzPred2" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# use the Stdz.DK.Biodiv object to generate a regression object named Stdz.Model1
Stdz.Model <- lm(---------- ~ ---------- + ---------- + ---------- + ---------- + ---------- + ---------- + ----------,
              data = ----------) # Call the object with the Scaled data

# use the summary() function to see the regression coefficients
----------(----------)
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="StdzPred2-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# use the Stdz.DK.Biodiv object to generate a regression object named Stdz.Model1
Stdz.Model <- lm(Log.S.AllGrp ~ NatDensBasMp + Dis2Cost + HII + Slope30mAgg + PrecSea + decimalLongitude + decimalLatitude,
              data = as.data.frame(Stdz.DK.Biodiv)) # Call the object with the Scaled data

# use the summary() function to see the regression coefficients
summary(Stdz.Model)
```

</div>


## 11. Multiple Linear Regressions - Summary of the Results

<div class="alert alert-info">
**Your task:**

Build a `data.frame` that contains the points in the table below. Do this by extracting information of the objects you have created beforehand, **NOT BY TYPING THE VALUES!**

+-------------+----------+----------------+--------------------------+-----------+---+---+
| Coefficient | Estimate | Standard Error | Standardised Coefficient | Tolerance | t |$p$|
+=============+==========+================+==========================+===========+===+===+
| Intercept   |          |                |                          |           |   |   |
+-------------+----------+----------------+--------------------------+-----------+---+---+
| LAT         |          |                |                          |           |   |   |
+-------------+----------+----------------+--------------------------+-----------+---+---+
| LONG        |          |                |                          |           |   |   |
+-------------+----------+----------------+--------------------------+-----------+---+---+
| MAP         |          |                |                          |           |   |   |
+-------------+----------+----------------+--------------------------+-----------+---+---+
| MAT         |          |                |                          |           |   |   |
+-------------+----------+----------------+--------------------------+-----------+---+---+
| JJAMAP      |          |                |                          |           |   |   |
+-------------+----------+----------------+--------------------------+-----------+---+---+
| DJFMAP      |          |                |                          |           |   |   |
+-------------+----------+----------------+--------------------------+-----------+---+---+

</div>

<div class="tutorial-exercise" data-label="SummTbl" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="20">

```text
# Generate three objects 
# 1. Store the summary of the regression carried out on original data (FullMod). name this sum.par.lm 
sum.par.lm <- ----------(----------)

#2.  Store the summary of the regression carried out on standardised data (Stdz.Model). name this sum.par.lm.scaled 
sum.par.lm.scaled <- ----------(----------)

#3. Create a vector called tolerances with the tolerance for each variable - use the function vif() for this
tolerances <- car::vif(Stdz.Model)

# Merge all the elements into a single data.frame
par.table <- cbind(----------, # un-standardised coefficients and SD
                   ----------, # standardized coefficients
                   c(NA, tolerances), # Tolerances
                   ----------) # un-standardised coefficients T value and P value
# Some housekeeping
colnames(par.table) <- c("Estimate", "Standard error", "Stdandard coefficient", "Tolerance", "t", "p")
round(par.table,3)
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="SummTbl-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Generate three objects 
# 1. Store the summary of the regression carried out on original data (FullMod). name this sum.par.lm 
sum.par.lm <- summary(FullMod)

#2.  Store the summary of the regression carried out on standardised data (Stdz.Model). name this sum.par.lm.scaled 
sum.par.lm.scaled <- summary(Stdz.Model)

#3. Create a vector called tolerances with the tolerance for each variable - use the function vif() for this
tolerances <- car::vif(Stdz.Model)

# Merge all the elements into a single data.frame
par.table <- cbind(sum.par.lm$coefficients[, 1:2], # un-standardised coefficients and SD
                   sum.par.lm.scaled$coefficients[, 1], # standardized coefficients
                   c(NA, tolerances), # Tolerances
                   sum.par.lm$coefficients[, 3:4]) # un-standardised coefficients T value and P value

colnames(par.table) <- c("Estimate", "Standard error", "Stdandard coefficient", "Tolerance", "t", "p")
round(par.table,3)
```

</div>

## 12. Multiple Linear Regressions - How good is my model

<div class="alert alert-info">
**Your task:**

Extract the model coefficient of determination and adjusted-coefficient of determination for both the modes. Do this for both the original model and the one run with standardised data.
  
</div>

<div class="tutorial-exercise" data-label="Rsqrd" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="10">

```text
# To extract the R^2 and adjusted-R^2 you need to use the `adj.r.squared` and `r.squared` subscripts on an object that builds the summary of an object created with the function `lm`
# r.squared for a regression carried out on original data 
----------(----------)$----------

# adj.r.squared for a regression carried out on original data 
----------(----------)$----------

# r.squared for a regression carried out on standardised data 
----------(----------)$----------

# adj.r.squared for a regression carried out on standardised data 
----------(----------)$----------
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="Rsqrd-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# To extract the R^2 and adjusted-R^2 you need to use the `adj.r.squared` and `r.squared` subscripts on an object that builds the summary of an object created with the function `lm`
# r.squared for a regression carried out on original data 
summary(FullMod)$r.squared

# adj.r.squared for a regression carried out on original data 
summary(FullMod)$adj.r.squared

# r.squared for a regression carried out on standardised data 
summary(Stdz.Model)$r.squared

# adj.r.squared for a regression carried out on standardised data 
summary(Stdz.Model)$adj.r.squared
```

</div>

```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="RsqrdQ1" class="tutorial-question panel-body">
<div id="RsqrdQ1-answer_container" class="shiny-html-output"></div>
<div id="RsqrdQ1-message_container" class="shiny-html-output"></div>
<div id="RsqrdQ1-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```

## 13. Multiple Linear Regressions - Minimum Adequate Model


<div class="alert alert-info">
**Your task:**

Run a forward selection procedure. For each process, extract the regression coefficients of the reduced model.
</div>

<div class="tutorial-exercise" data-label="MAM1" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="10">

```text
# Using the step() function, make a forward stepwise selection on the model with standardised regression coefficients. Save this into a object call FrwdSel
FrwdSel <- ----------(----------, # an object representing a model to simplify
                direction = "----------") # forward or backward

#extract the regression coefficients using the function coef()
----------(----------)

#Get the adjusted R-squared
----------(----------)$----------
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>

<div class="tutorial-exercise-support" data-label="MAM1-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Using the step() function, make a forward stepwise selection on the model with standardised regression coefficients. Save this into a object call FrwdSel
FrwdSel <- step(Stdz.Model, # an object representing a model to simplify
                direction = "forward") # forward or backward
#extract the regression coefficients using the function coef()
coef(FrwdSel)
#Get the adjusted R-squared
summary(FrwdSel)$adj.r.squared
```

</div>

<div class="alert alert-info">
**Your task:**

Run a backward selection procedure. For each process, extract the regression coefficients of the reduced model.
</div>

<div class="tutorial-exercise" data-label="MAM2" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="10">

```text
# Using the step() function, make a Backwards stepwise selection on the model with standardised regression coefficients. Save this into a object call BackSel
BackSel <- ----------(----------, # an object representing a model to simplify
                direction = "----------") # forward or backward

#extract the regression coefficients using the function coef()
----------(----------)

#Get the adjusted R-squared
----------(----------)$----------
```

<script type="application/json" data-ui-opts="1">{"engine":"r","has_checker":false,"caption":"<span data-i18n=\"text.enginecap\" data-i18n-opts=\"{&quot;engine&quot;:&quot;R&quot;}\">R Code<\/span>"}</script></div>
<div class="tutorial-exercise-support" data-label="MAM2-solution" data-completion="1" data-diagnostics="1" data-startover="1" data-lines="0">

```text
# Using the step() function, make a Backwards stepwise selection on the model with standardised regression coefficients. Save this into a object call BackSel
BackSel <- step(Stdz.Model, # an object representing a model to simplify
                direction = "backward",
                na.rm=T) # forward or backward
#extract the regression coefficients using the function coef()
coef(BackSel)
#Get the adjusted R-squared
summary(BackSel)$adj.r.squared
```

</div>


```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="MAMQ1" class="tutorial-question panel-body">
<div id="MAMQ1-answer_container" class="shiny-html-output"></div>
<div id="MAMQ1-message_container" class="shiny-html-output"></div>
<div id="MAMQ1-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```

```{=html}
<div class="panel panel-default tutorial-question-container">
<div data-label="MAMQ2" class="tutorial-question panel-body">
<div id="MAMQ2-answer_container" class="shiny-html-output"></div>
<div id="MAMQ2-message_container" class="shiny-html-output"></div>
<div id="MAMQ2-action_button_container" class="shiny-html-output"></div>
<script>if (Tutorial.triggerMathJax) Tutorial.triggerMathJax()</script>
</div>
</div>
```
preserve290b29455a56f413
preserve473a8476e6a7026d
preservef019b45861b7578a
preserve2ed1c7cc4477c5aa
preserve946a495e6d712c07
preserve14bb99aa48ec4321
preserve69eb434f33568093
preserve4b300f879522bd41
preserve48c8c3563d883397
preservef6b73202e11cd03f
preserve81369397cd1269da
preserve6eb5612db79d8dd5
preserve753a15d84e3e38ea
preserve372428dcecaf7888
preserve6016c644c9146cff
preservea72f9c8a0f816f00
preserve413bcb5a80206d12
preserveaa7e581d779a94bb
preserve7aa28939dd9db0cb
preserve1f1775646f4a060a
preserve1aa483709763a299
preserve2078d4a3b1bb8ca7
preservebebac8e1dc368ea1
preserveddbfe29bedbd1005
preserve51d7e7c030a0e7b7
preserve8d4e59451ca5b545
preserve4d77580ace3f4625
preserve840298f814cb74a1
preserve8b39e7233194649c
preserve72a0a1dfea44a72d
preserveb4f153d095314a30
preserve49ca51dc9bc93d29
preserve5a6e7898d6288409
preserve8d7ec7008ed072fa
preservecf9e58535ec56c83
preserved8a09b80259363a5
preserve705bcd1fe6bcf7e0
preserve76744a6470a48eed
preservebb58afd9eb533c74
preservea186d63b664144fb
preservea8f02b5bff5c1ab8
preserve4727381050f0ca3e
preservea4d1496652ed29ad
preservebdb5b4156f1f8526
preservec85ae91c5cbaff54
preserve4ff72cd00c87111a
preserve5bd6f58276e09e6d
preserve2c7213ca70dbcb9d
preserveda61b588a0facc38
preserveec9282460603b65c
preserve921a78882a91f53c

<!--html_preserve-->
<script type="application/shiny-prerendered" data-context="dependencies">
{"type":"list","attributes":{},"value":[{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["header-attrs"]},{"type":"character","attributes":{},"value":["2.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["rmd/h/pandoc"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["header-attrs.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["rmarkdown"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["2.14"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["jquery"]},{"type":"character","attributes":{},"value":["3.6.0"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/3.6.0"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["jquery-3.6.0.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["jquerylib"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.1.4"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["bootstrap"]},{"type":"character","attributes":{},"value":["3.3.5"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["rmd/h/bootstrap"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["viewport"]}},"value":[{"type":"character","attributes":{},"value":["width=device-width, initial-scale=1"]}]},{"type":"character","attributes":{},"value":["js/bootstrap.min.js","shim/html5shiv.min.js","shim/respond.min.js"]},{"type":"character","attributes":{},"value":["css/cerulean.min.css"]},{"type":"character","attributes":{},"value":["<style>h1 {font-size: 34px;}\n       h1.title {font-size: 38px;}\n       h2 {font-size: 30px;}\n       h3 {font-size: 24px;}\n       h4 {font-size: 18px;}\n       h5 {font-size: 16px;}\n       h6 {font-size: 12px;}\n       code {color: inherit; background-color: rgba(0, 0, 0, 0.04);}\n       pre:not([class]) { background-color: white }<\/style>"]},{"type":"NULL"},{"type":"character","attributes":{},"value":["rmarkdown"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["2.14"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["pagedtable"]},{"type":"character","attributes":{},"value":["1.1"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["rmd/h/pagedtable-1.1"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["js/pagedtable.js"]},{"type":"character","attributes":{},"value":["css/pagedtable.css"]},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["rmarkdown"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["2.14"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["highlightjs"]},{"type":"character","attributes":{},"value":["9.12.0"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["rmd/h/highlightjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["highlight.js"]},{"type":"character","attributes":{},"value":["textmate.css"]},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["rmarkdown"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["2.14"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["tutorial"]},{"type":"character","attributes":{},"value":["0.10.5.9000"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/tutorial"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["tutorial.js"]},{"type":"character","attributes":{},"value":["tutorial.css"]},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["i18n"]},{"type":"character","attributes":{},"value":["21.6.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/i18n"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["i18next.min.js","tutorial-i18n-init.js"]},{"type":"NULL"},{"type":"character","attributes":{},"value":["<script id=\"i18n-cstm-trns\" type=\"application/json\">{\"language\":\"en\",\"resources\":{\"en\":{\"translation\":{\"button\":{\"runcode\":\"Run Code\",\"runcodetitle\":\"$t(button.runcode) ({{kbd}})\",\"hint\":\"Hint\",\"hint_plural\":\"Hints\",\"hinttitle\":\"$t(button.hint)\",\"hintnext\":\"Next Hint\",\"hintprev\":\"Previous Hint\",\"solution\":\"Solution\",\"solutiontitle\":\"$t(button.solution)\",\"copyclipboard\":\"Copy to Clipboard\",\"startover\":\"Start Over\",\"startovertitle\":\"$t(button.startover)\",\"continue\":\"Continue\",\"submitanswer\":\"Submit Answer\",\"submitanswertitle\":\"$t(button.submitanswer)\",\"previoustopic\":\"Previous Topic\",\"nexttopic\":\"Next Topic\",\"questionsubmit\":\"$t(button.submitanswer)\",\"questiontryagain\":\"Try Again\"},\"text\":{\"startover\":\"Start Over\",\"areyousure\":\"Are you sure you want to start over? (all exercise progress will be reset)\",\"youmustcomplete\":\"You must complete the\",\"exercise\":\"exercise\",\"exercise_plural\":\"exercises\",\"inthissection\":\"in this section before continuing.\",\"code\":\"Code\",\"enginecap\":\"{{engine}} $t(text.code)\",\"quiz\":\"Quiz\",\"blank\":\"blank\",\"blank_plural\":\"blanks\",\"exercisecontainsblank\":\"This exercise contains {{count}} $t(text.blank).\",\"pleasereplaceblank\":\"Please replace {{blank}} with valid code.\",\"unparsable\":\"It looks like this might not be valid R code. R cannot determine how to turn your text into a complete command. You may have forgotten to fill in a blank, to remove an underscore, to include a comma between arguments, or to close an opening <code>&quot;<\\/code>, <code>'<\\/code>, <code>(<\\/code> or <code>{<\\/code> with a matching <code>&quot;<\\/code>, <code>'<\\/code>, <code>)<\\/code> or <code>}<\\/code>.\\n\",\"unparsablequotes\":\"<p>It looks like your R code contains specially formatted quotation marks or &quot;curly&quot; quotes (<code>{{character}}<\\/code>) around character strings, making your code invalid. R requires character values to be contained in straight quotation marks (<code>&quot;<\\/code> or <code>'<\\/code>).<\\/p> {{code}} <p>Don't worry, this is a common source of errors when you copy code from another app that applies its own formatting to text. You can try replacing the code on that line with the following. There may be other places that need to be fixed, too.<\\/p> {{suggestion}}\\n\",\"unparsableunicode\":\"<p>It looks like your R code contains an unexpected special character (<code>{{character}}<\\/code>) that makes your code invalid.<\\/p> {{code}} <p>Sometimes your code may contain a special character that looks like a regular character, especially if you copy and paste the code from another app. Try deleting the special character from your code and retyping it manually.<\\/p>\\n\",\"unparsableunicodesuggestion\":\"<p>It looks like your R code contains an unexpected special character (<code>{{character}}<\\/code>) that makes your code invalid.<\\/p> {{code}} <p>Sometimes your code may contain a special character that looks like a regular character, especially if you copy and paste the code from another app. You can try replacing the code on that line with the following. There may be other places that need to be fixed, too.<\\/p> {{suggestion}}\\n\",\"and\":\"and\",\"or\":\"or\",\"listcomma\":\", \",\"oxfordcomma\":\",\"}}},\"fr\":{\"translation\":{\"button\":{\"runcode\":\"Lancer le Code\",\"runcodetitle\":\"$t(button.runcode) ({{kbd}})\",\"hint\":\"Indication\",\"hint_plural\":\"Indications\",\"hinttitle\":\"$t(button.hint)\",\"hintnext\":\"Indication Suivante\",\"hintprev\":\"Indication Précédente\",\"solution\":\"Solution\",\"solutiontitle\":\"$t(button.solution)\",\"copyclipboard\":\"Copier dans le Presse-papier\",\"startover\":\"Recommencer\",\"startovertitle\":\"$t(button.startover)\",\"continue\":\"Continuer\",\"submitanswer\":\"Soumettre\",\"submitanswertitle\":\"$t(button.submitanswer)\",\"previoustopic\":\"Chapitre Précédent\",\"nexttopic\":\"Chapitre Suivant\",\"questionsubmit\":\"$t(button.submitanswer)\",\"questiontryagain\":\"Réessayer\"},\"text\":{\"startover\":\"Recommencer\",\"areyousure\":\"Êtes-vous certains de vouloir recommencer? (La progression sera remise à zéro)\",\"youmustcomplete\":\"Vous devez d'abord compléter\",\"exercise\":\"l'exercice\",\"exercise_plural\":\"des exercices\",\"inthissection\":\"de cette section avec de continuer.\",\"code\":\"Code\",\"enginecap\":\"$t(text.code) {{engine}}\",\"quiz\":\"Quiz\",\"and\":\"et\",\"or\":\"ou\",\"oxfordcomma\":\"\"}}},\"es\":{\"translation\":{\"button\":{\"runcode\":\"Ejecutar código\",\"runcodetitle\":\"$t(button.runcode) ({{kbd}})\",\"hint\":\"Pista\",\"hint_plural\":\"Pistas\",\"hinttitle\":\"$t(button.hint)\",\"hintnext\":\"Siguiente pista\",\"hintprev\":\"Pista anterior\",\"solution\":\"Solución\",\"solutiontitle\":\"$t(button.solution)\",\"copyclipboard\":\"Copiar al portapapeles\",\"startover\":\"Reiniciar\",\"startovertitle\":\"$t(button.startover)\",\"continue\":\"Continuar\",\"submitanswer\":\"Enviar respuesta\",\"submitanswertitle\":\"$t(button.submitanswer)\",\"previoustopic\":\"Tema anterior\",\"nexttopic\":\"Tema siguiente\",\"questionsubmit\":\"$t(button.submitanswer)\",\"questiontryagain\":\"Volver a intentar\"},\"text\":{\"startover\":\"Reiniciar\",\"areyousure\":\"¿De verdad quieres empezar de nuevo? (todo el progreso del ejercicio se perderá)\",\"youmustcomplete\":\"Debes completar\",\"exercise\":\"el ejercicio\",\"exercise_plural\":\"los ejercicios\",\"inthissection\":\"en esta sección antes de continuar.\",\"code\":\"Código\",\"enginecap\":\"$t(text.code) {{engine}}\",\"quiz\":\"Cuestionario\",\"and\":\"y\",\"or\":\"o\",\"oxfordcomma\":\"\"}}},\"pt\":{\"translation\":{\"button\":{\"runcode\":\"Executar código\",\"runcodetitle\":\"$t(button.runcode) ({{kbd}})\",\"hint\":\"Dica\",\"hint_plural\":\"Dicas\",\"hinttitle\":\"$t(button.hint)\",\"hintnext\":\"Próxima dica\",\"hintprev\":\"Dica anterior\",\"solution\":\"Solução\",\"solutiontitle\":\"$t(button.solution)\",\"copyclipboard\":\"Copiar para a área de transferência\",\"startover\":\"Reiniciar\",\"startovertitle\":\"$t(button.startover)\",\"continue\":\"Continuar\",\"submitanswer\":\"Enviar resposta\",\"submitanswertitle\":\"$t(button.submitanswer)\",\"previoustopic\":\"Tópico anterior\",\"nexttopic\":\"Próximo tópico\",\"questionsubmit\":\"$t(button.submitanswer)\",\"questiontryagain\":\"Tentar novamente\"},\"text\":{\"startover\":\"Reiniciar\",\"areyousure\":\"Tem certeza que deseja começar novamente? (todo o progresso feito será perdido)\",\"youmustcomplete\":\"Você deve completar\",\"exercise\":\"o exercício\",\"exercise_plural\":\"os exercícios\",\"inthissection\":\"nesta seção antes de continuar.\",\"code\":\"Código\",\"enginecap\":\"$t(text.code) {{engine}}\",\"quiz\":\"Quiz\",\"and\":\"e\",\"or\":\"ou\",\"oxfordcomma\":\"\"}}},\"tr\":{\"translation\":{\"button\":{\"runcode\":\"Çalıştırma Kodu\",\"runcodetitle\":\"$t(button.runcode) ({{kbd}})\",\"hint\":\"Ipucu\",\"hint_plural\":\"İpuçları\",\"hinttitle\":\"$t(button.hint)\",\"hintnext\":\"Sonraki İpucu\",\"hintprev\":\"Önceki İpucu\",\"solution\":\"Çözüm\",\"solutiontitle\":\"$t(button.solution)\",\"copyclipboard\":\"Pano'ya Kopyala\",\"startover\":\"Baştan Başlamak\",\"startovertitle\":\"$t(button.startover)\",\"continue\":\"Devam et\",\"submitanswer\":\"Cevabı onayla\",\"submitanswertitle\":\"$t(button.submitanswer)\",\"previoustopic\":\"Önceki Konu\",\"nexttopic\":\"Sonraki Konu\",\"questionsubmit\":\"$t(button.submitanswer)\",\"questiontryagain\":\"Tekrar Deneyin\"},\"text\":{\"startover\":\"Baştan Başlamak\",\"areyousure\":\"Baştan başlamak istediğinizden emin misiniz? (tüm egzersiz ilerlemesi kaybolacak)\",\"youmustcomplete\":\"Tamamlamalısın\",\"exercise\":\"egzersiz\",\"exercise_plural\":\"egzersizler\",\"inthissection\":\"devam etmeden önce bu bölümde\",\"code\":\"Kod\",\"enginecap\":\"$t(text.code) {{engine}}\",\"quiz\":\"Sınav\",\"oxfordcomma\":\"\"}}},\"emo\":{\"translation\":{\"button\":{\"runcode\":\"🏃\",\"runcodetitle\":\"$t(button.runcode) ({{kbd}})\",\"hint\":\"💡\",\"hint_plural\":\"$t(button.hint)\",\"hinttitle\":\"$t(button.hint)\",\"solution\":\"🎯\",\"solutiontitle\":\"$t(button.solution)\",\"copyclipboard\":\"📋\",\"startover\":\"⏮\",\"startovertitle\":\"Start Over\",\"continue\":\"✅\",\"submitanswer\":\"🆗\",\"submitanswertitle\":\"Submit Answer\",\"previoustopic\":\"⬅\",\"nexttopic\":\"➡\",\"questionsubmit\":\"$t(button.submitanswer)\",\"questiontryagain\":\"🔁\"},\"text\":{\"startover\":\"⏮\",\"areyousure\":\"🤔\",\"youmustcomplete\":\"⚠️ 👉 🧑‍💻\",\"exercise\":\"\",\"exercise_plural\":\"\",\"inthissection\":\"\",\"code\":\"💻\",\"enginecap\":\"$t(text.code) {{engine}}\",\"oxfordcomma\":\"\"}}},\"eu\":{\"translation\":{\"button\":{\"runcode\":\"Kodea egikaritu\",\"runcodetitle\":\"$t(button.runcode) ({{kbd}})\",\"hint\":\"Laguntza\",\"hint_plural\":\"Laguntza\",\"hinttitle\":\"$t(button.hint)\",\"hintnext\":\"Aurreko laguntza\",\"hintprev\":\"Hurrengo laguntza\",\"solution\":\"Ebazpena\",\"solutiontitle\":\"$t(button.solution)\",\"copyclipboard\":\"Arbelean kopiatu\",\"startover\":\"Berrabiarazi\",\"startovertitle\":\"$t(button.startover)\",\"continue\":\"Jarraitu\",\"submitanswer\":\"Erantzuna bidali\",\"submitanswertitle\":\"$t(button.submitanswer)\",\"previoustopic\":\"Aurreko atala\",\"nexttopic\":\"Hurrengo atala\",\"questionsubmit\":\"$t(button.submitanswer)\",\"questiontryagain\":\"Berriro saiatu\"},\"text\":{\"startover\":\"Berrabiarazi\",\"areyousure\":\"Berriro hasi nahi duzu? (egindako lana galdu egingo da)\",\"youmustcomplete\":\"Aurrera egin baino lehen atal honetako\",\"exercise\":\"ariketa egin behar duzu.\",\"exercise_plural\":\"ariketak egin behar dituzu.\",\"inthissection\":\"\",\"code\":\"Kodea\",\"enginecap\":\"$t(text.code) {{engine}}\",\"quiz\":\"Galdetegia\",\"oxfordcomma\":\"\"}}},\"de\":{\"translation\":{\"button\":{\"runcode\":\"Code ausführen\",\"runcodetitle\":\"$t(button.runcode) ({{kbd}})\",\"hint\":\"Tipp\",\"hint_plural\":\"Tipps\",\"hinttitle\":\"$t(button.hint)\",\"hintnext\":\"Nächster Tipp\",\"hintprev\":\"Vorheriger Tipp\",\"solution\":\"Lösung\",\"solutiontitle\":\"$t(button.solution)\",\"copyclipboard\":\"In die Zwischenablage kopieren\",\"startover\":\"Neustart\",\"startovertitle\":\"$t(button.startover)\",\"continue\":\"Weiter\",\"submitanswer\":\"Antwort einreichen\",\"submitanswertitle\":\"$t(button.submitanswer)\",\"previoustopic\":\"Vorheriges Kapitel\",\"nexttopic\":\"Nächstes Kapitel\",\"questionsubmit\":\"$t(button.submitanswer)\",\"questiontryagain\":\"Nochmal versuchen\"},\"text\":{\"startover\":\"Neustart\",\"areyousure\":\"Bist du sicher, dass du neustarten willst? (der gesamte Lernfortschritt wird gelöscht)\",\"youmustcomplete\":\"Vervollstädinge\",\"exercise\":\"die Übung\",\"exercise_plural\":\"die Übungen\",\"inthissection\":\"in diesem Kapitel, bevor du fortfährst.\",\"code\":\"Code\",\"enginecap\":\"$t(text.code) {{engine}}\",\"quiz\":\"Quiz\",\"blank\":\"Lücke\",\"blank_plural\":\"Lücken\",\"pleasereplaceblank\":\"Bitte ersetze {{blank}} mit gültigem Code.\",\"unparsable\":\"Dies scheint kein gültiger R Code zu sein. R kann deinen Text nicht in einen gültigen Befehl übersetzen. Du hast vielleicht vergessen, die Lücke zu füllen, einen Unterstrich zu entfernen, ein Komma zwischen Argumente zu setzen oder ein eröffnendes <code>&quot;<\\/code>, <code>'<\\/code>, <code>(<\\/code> oder <code>{<\\/code> mit einem zugehörigen <code>&quot;<\\/code>, <code>'<\\/code>, <code>)<\\/code> oder <code>}<\\/code> zu schließen.\\n\",\"and\":\"und\",\"or\":\"oder\",\"listcomma\":\", \",\"oxfordcomma\":\",\"}}},\"ko\":{\"translation\":{\"button\":{\"runcode\":\"코드 실행\",\"runcodetitle\":\"$t(button.runcode) ({{kbd}})\",\"hint\":\"힌트\",\"hint_plural\":\"힌트들\",\"hinttitle\":\"$t(button.hint)\",\"hintnext\":\"다음 힌트\",\"hintprev\":\"이전 힌트\",\"solution\":\"솔루션\",\"solutiontitle\":\"$t(button.solution)\",\"copyclipboard\":\"클립보드에 복사\",\"startover\":\"재학습\",\"startovertitle\":\"$t(button.startover)\",\"continue\":\"다음 학습으로\",\"submitanswer\":\"정답 제출\",\"submitanswertitle\":\"$t(button.submitanswer)\",\"previoustopic\":\"이전 토픽\",\"nexttopic\":\"다음 토픽\",\"questionsubmit\":\"$t(button.submitanswer)\",\"questiontryagain\":\"재시도\"},\"text\":{\"startover\":\"재학습\",\"areyousure\":\"다시 시작 하시겠습니까? (모든 예제의 진행 정보가 재설정됩니다)\",\"youmustcomplete\":\"당신은 완료해야 합니다\",\"exercise\":\"연습문제\",\"exercise_plural\":\"연습문제들\",\"inthissection\":\"이 섹션을 실행하기 전에\",\"code\":\"코드\",\"enginecap\":\"$t(text.code) {{engine}}\",\"quiz\":\"퀴즈\",\"blank\":\"공백\",\"blank_plural\":\"공백들\",\"exercisecontainsblank\":\"이 연습문제에는 {{count}}개의 $t(text.blank)이 포함되어 있습니다.\",\"pleasereplaceblank\":\"{{blank}}를 유효한 코드로 바꾸십시오.\",\"unparsable\":\"이것은 유효한 R 코드가 아닐 수 있습니다. R은 텍스트를 완전한 명령으로 변환하는 방법을 결정할 수 없습니다. 당신은 공백이나 밑줄을 대체하여 채우기, 인수를 컴마로 구분하기, 또는 <code>&quot;<\\/code>, <code>'<\\/code>, <code>(<\\/code> , <code>{<\\/code>로 시작하는 구문을 닫는 <code>&quot;<\\/code>, <code>'<\\/code>, <code>)<\\/code>, <code>}<\\/code>을 잊었을 수도 있습니다.\\n\",\"and\":\"그리고\",\"or\":\"혹은\",\"listcomma\":\", \",\"oxfordcomma\":\"\"}}},\"zh\":{\"translation\":{\"button\":{\"runcode\":\"运行代码\",\"runcodetitle\":\"$t(button.runcode) ({{kbd}})\",\"hint\":\"提示\",\"hint_plural\":\"提示\",\"hinttitle\":\"$t(button.hint)\",\"hintnext\":\"下一个提示\",\"hintprev\":\"上一个提示\",\"solution\":\"答案\",\"solutiontitle\":\"$t(button.solution)\",\"copyclipboard\":\"复制到剪切板\",\"startover\":\"重新开始\",\"startovertitle\":\"$t(button.startover)\",\"continue\":\"继续\",\"submitanswer\":\"提交答案\",\"submitanswertitle\":\"$t(button.submitanswer)\",\"previoustopic\":\"上一专题\",\"nexttopic\":\"下一专题\",\"questionsubmit\":\"$t(button.submitanswer)\",\"questiontryagain\":\"再试一次\"},\"text\":{\"startover\":\"重置\",\"areyousure\":\"你确定要重新开始吗? (所有当前进度将被重置)\",\"youmustcomplete\":\"你必须完成\",\"exercise\":\"练习\",\"exercise_plural\":\"练习\",\"inthissection\":\"在进行本节之前\",\"code\":\"代码\",\"enginecap\":\"$t(text.code) {{engine}}\",\"quiz\":\"测试\",\"blank\":\"空\",\"blank_plural\":\"空\",\"exercisecontainsblank\":\"本练习包含{{count}}个$t(text.blank)\",\"pleasereplaceblank\":\"请在{{blank}}内填写恰当的代码\",\"unparsable\":\"这似乎不是有效的R代码。 R不知道如何将您的文本转换为完整的命令。 您是否忘了填空，忘了删除下划线，忘了在参数之间包含逗号，或者是忘了用<code>&quot;<\\/code>, <code>'<\\/code>, <code>)<\\/code>,<code>}<\\/code>来封闭<code>&quot;<\\/code>, <code>'<\\/code>, <code>(<\\/code>。 or <code>{<\\/code>。\\n\",\"unparsablequotes\":\"<p>您的R代码中似乎含有特殊格式的引号，或者弯引号(<code>{{character}}<\\/code>) 在字符串前后，在R中字符串应该被直引号(<code>&quot;<\\/code> 或者 <code>'<\\/code>)包裹。<\\/p> {{code}} <p>别担心，该错误经常在复制粘贴包含格式的代码时遇到， 您可以尝试将该行中的代码替换为以下代码，也许还有其他地方需要修改。<\\/p> {{suggestion}}\\n\",\"unparsableunicode\":\"<p>您的代码中似乎包含有异常字符(<code>{{character}}<\\/code>),导致代码无效。<\\/p> {{code}} <p>有时候你的代码可能含有看似正常字符的特殊字符，特别是当你复制粘贴其他来源代码的时候。 请试着删除这些特殊字符,重新输入<\\/p>\\n\",\"unparsableunicodesuggestion\":\"<p>您的代码中似乎包含有异常字符(<code>{{character}}<\\/code>),导致代码无效。<\\/p> {{code}} <p>有时候你的代码可能含有看似正常字符的特殊字符，特别是当你复制粘贴其他来源代码的时候。 请试着删除这些特殊字符,重新输入<\\/p>\\n\",\"and\":\"且\",\"or\":\"或\",\"listcomma\":\",\",\"oxfordcomma\":\",\"}}},\"pl\":{\"translation\":{\"button\":{\"runcode\":\"Uruchom kod\",\"runcodetitle\":\"$t(button.runcode) ({{kbd}})\",\"hint\":\"Podpowiedź\",\"hint_plural\":\"Podpowiedzi\",\"hinttitle\":\"$t(button.hint)\",\"hintnext\":\"Następna podpowiedź\",\"hintprev\":\"Poprzednia podpowiedź\",\"solution\":\"Rozwiązanie\",\"solutiontitle\":\"$t(button.solution)\",\"copyclipboard\":\"Kopiuj do schowka\",\"startover\":\"Zacznij od początku\",\"startovertitle\":\"$t(button.startover)\",\"continue\":\"Kontynuuj\",\"submitanswer\":\"Wyślij\",\"submitanswertitle\":\"$t(button.submitanswer)\",\"previoustopic\":\"Poprzednia sekcja\",\"nexttopic\":\"Następna sekcja\",\"questionsubmit\":\"$t(button.submitanswer)\",\"questiontryagain\":\"Spróbuj ponownie\"},\"text\":{\"startover\":\"Zacznij od początku\",\"areyousure\":\"Czy na pewno chcesz zacząć od początku? (cały postęp w zadaniu zostanie utracony)\",\"youmustcomplete\":\"Musisz ukończyć\",\"exercise\":\"ćwiczenie\",\"exercise_plural\":\"ćwiczenia\",\"inthissection\":\"w tej sekcji przed kontynuowaniem\",\"code\":\"Kod\",\"enginecap\":\"$t(text.code) {{engine}}\",\"quiz\":\"Quiz\",\"blank\":\"luka\",\"blank_plural\":\"luk(i)\",\"exercisecontainsblank\":\"To ćwiczenie zawiera {{count}} $t(text.blank).\",\"pleasereplaceblank\":\"Proszę uzupełnić {{blank}} prawidłowym kodem.\",\"unparsable\":\"Wygląda na to, że może to nie być prawidłowy kod R. R nie jest w stanie przetworzyć Twojego tekstu na polecenie. Mogłeś(-aś) zapomnieć wypełnić luki, usunąć podkreślnik, umieścić przecinka między argumentami, lub zamknąć znak <code>&quot;<\\/code>, <code>'<\\/code>, <code>(<\\/code> lub <code>{<\\/code> odpowiadającym <code>&quot;<\\/code>, <code>'<\\/code>, <code>)<\\/code> lub <code>}<\\/code>.\\n\",\"unparsablequotes\":\"<p>Wygląda na to, że Twój kod zawiera szczególnie sformatowane cudzysłowy lub cudzysłowy typograficzne (<code>{{character}}<\\/code>) przy ciągach znaków, co sprawia, że kod jest niepoprawny. R wymaga cudzysłowów prostych (<code>&quot;<\\/code> albo <code>'<\\/code>).<\\/p> {{code}} <p>Nie martw się, to powszechne źródło błędów, gdy kopiuje się kod z innego programu, który sam formatuje teskt. Możesz spróbować zastąpić swój kod następującym kodem. Mogą być też inne miejsca, które wymagają poprawienia.<\\/p> {{suggestion}}\\n\",\"unparsableunicode\":\"<p>Wygląda na to, że Twój kod zawiera niespodziewany znak specjalny (<code>{{character}}<\\/code>), co sprawia, że kod jest niepoprawny.<\\/p> {{code}} <p>Czasami Twój kod może zawierać znak specjalny, który wygląda jak zwykły znak, zwłaszcza jeśli kopiujesz kod z innego programu. Spróbuj usunąć znak specjalny i wpisać do ponownie ręcznie.<\\/p>\\n\",\"unparsableunicodesuggestion\":\"<p>Wygląda na to, że Twój kod zawiera niespodziewany znak specjalny (<code>{{character}}<\\/code>), co sprawia, że kod jest niepoprawny.<\\/p> {{code}} <p>Czasami Twój kod może zawierać znak specjalny, który wygląda jak zwykły znak, zwłaszcza jeśli kopiujesz kod z innego programu. Możesz spróbować zastąpić swój kod następującym kodem. Mogą być też inne miejsca, które wymagają poprawienia.<\\/p> {{suggestion}}\\n\",\"and\":\"i\",\"or\":\"lub\",\"listcomma\":\", \",\"oxfordcomma\":\"\"}}}}}<\/script>"]},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["tutorial-format"]},{"type":"character","attributes":{},"value":["0.10.5.9000"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["rmarkdown/templates/tutorial/resources"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["tutorial-format.js"]},{"type":"character","attributes":{},"value":["tutorial-format.css","rstudio-theme.css"]},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["jquery"]},{"type":"character","attributes":{},"value":["3.6.0"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/3.6.0"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["jquery-3.6.0.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["jquerylib"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.1.4"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["navigation"]},{"type":"character","attributes":{},"value":["1.1"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["rmd/h/navigation-1.1"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["tabsets.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["rmarkdown"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["2.14"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["highlightjs"]},{"type":"character","attributes":{},"value":["9.12.0"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["rmd/h/highlightjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["highlight.js"]},{"type":"character","attributes":{},"value":["default.css"]},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["rmarkdown"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["2.14"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["jquery"]},{"type":"character","attributes":{},"value":["3.6.0"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/3.6.0"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["jquery-3.6.0.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["jquerylib"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.1.4"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["font-awesome"]},{"type":"character","attributes":{},"value":["5.1.0"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["rmd/h/fontawesome"]}]},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["css/all.css","css/v4-shims.css"]},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["rmarkdown"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["2.14"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["bootbox"]},{"type":"character","attributes":{},"value":["5.5.2"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/bootbox"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["bootbox.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["idb-keyvalue"]},{"type":"character","attributes":{},"value":["3.2.0"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/idb-keyval"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["idb-keyval-iife-compat.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[false]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["tutorial"]},{"type":"character","attributes":{},"value":["0.10.5.9000"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/tutorial"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["tutorial.js"]},{"type":"character","attributes":{},"value":["tutorial.css"]},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["ace"]},{"type":"character","attributes":{},"value":["1.4.14"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/ace"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["ace.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["name","version","src","meta","script","stylesheet","head","attachment","package","all_files","pkgVersion"]},"class":{"type":"character","attributes":{},"value":["html_dependency"]}},"value":[{"type":"character","attributes":{},"value":["clipboardjs"]},{"type":"character","attributes":{},"value":["2.0.10"]},{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["file"]}},"value":[{"type":"character","attributes":{},"value":["lib/clipboardjs"]}]},{"type":"NULL"},{"type":"character","attributes":{},"value":["clipboard.min.js"]},{"type":"NULL"},{"type":"NULL"},{"type":"NULL"},{"type":"character","attributes":{},"value":["learnr"]},{"type":"logical","attributes":{},"value":[true]},{"type":"character","attributes":{},"value":["0.10.5.9000"]}]}]}
</script>
<!--/html_preserve-->
<!--html_preserve-->
<script type="application/shiny-prerendered" data-context="execution_dependencies">
{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["packages"]}},"value":[{"type":"list","attributes":{"names":{"type":"character","attributes":{},"value":["packages","version"]},"class":{"type":"character","attributes":{},"value":["data.frame"]},"row.names":{"type":"integer","attributes":{},"value":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46]}},"value":[{"type":"character","attributes":{},"value":["backports","base","bslib","cachem","checkmate","cli","compiler","curl","datasets","digest","ellipsis","evaluate","fastmap","graphics","grDevices","htmltools","htmlwidgets","httpuv","jquerylib","jsonlite","knitr","later","learnr","lifecycle","magrittr","markdown","methods","mime","promises","R6","Rcpp","rlang","rmarkdown","rprojroot","rstudioapi","sass","shiny","stats","stringi","stringr","tools","utils","withr","xfun","xtable","yaml"]},{"type":"character","attributes":{},"value":["1.4.1","4.2.1","0.4.0","1.0.6","2.1.0","3.3.0","4.2.1","4.3.2","4.2.1","0.6.29","0.3.2","0.15","1.1.0","4.2.1","4.2.1","0.5.3","1.5.4","1.6.5","0.1.4","1.8.0","1.39","1.3.0","0.10.5.9000","1.0.1","2.0.3","1.1","4.2.1","0.12","1.2.0.1","2.5.1","1.0.9","1.0.4","2.14","2.0.3","0.13","0.4.2","1.7.2","4.2.1","1.7.8","1.4.0","4.2.1","4.2.1","2.5.0","0.31","1.8-4","2.3.5"]}]}]}
</script>
<!--/html_preserve-->
