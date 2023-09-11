What is in here:
Paruelo & Lauenroth (1996) analyzed the geographic distribution and the effects of climate variables on the relative abundance of a number of plant functional types (PFTs) including shrubs, forbs, succulents (e.g. cacti), C3 grasses and C4 grasses. The latter two PFTs represent grasses that utilize the C from the atmosphere differently in photosynthesis and are expected to have different responses to CO2 and climate change. They used data from 73 sites across temperate central North America and calculated the relative abundance of a series of plant functional types (PFTs) including C3 grasses and C4 grasses. These two PFTs represent grasses that utilize the C from the atmosphere differently in photosynthesis and are expected to have different responses to CO2 and climate change.

They analyzed the geographic distribution and the effects of climate variables on the relative abundance of C3 and C4 grasses. **Here you will only focus on C3 plants** (`LC3`), which relative abundance was positively skewed and transformed using a log~10~(x+0.1) transformation (log~10~(C3+0.1)).

Predictor variables used to describe each evaluated site included the latitude in centesimal degrees (`LAT`), the longitude in centesimal degrees (`LONG`), the mean annual precipitation in mm (`MAP`), the mean annual temperature in C (`MAT`), the proportion of `MAP` that fell in June, July and August (`JJAMAP`) and the proportion of `MAP` that fell in December, January and February (`DJFMAP`).

Format of the files:
	* pareulo.csv (comma delimeted ascii text file)

File contents:
	* C3 - relative abundance of C3 grasses
	* C4 - relative abundance of C4 grasses
	* MAP - mean annual precipitation (mm)
	* MAT - mean annual temperature (oC)
	* JJAMAP - proportion of MAP that fell in June, July and August
	* DJFMAP - proportion of MAP that fell in December, January and February
	* LONG - longitude in centesimal degrees
	* LAT - latitude in centesimal degrees
	* LC3 - log10 transformation of C3
	* LC4 - log10 transformation of C4
	* CLONG - centered LONG
	* CLAT - centered LAT
	* RESID1 - residuals from linear regression of LC3 against CLAT + CLONG + CLAT*CLONG
	* PREDICT1 - predicted LC3 from linear regression of LC3 against CLAT + CLONG + CLAT*CLONG
----------------------------------------------------------------------------------