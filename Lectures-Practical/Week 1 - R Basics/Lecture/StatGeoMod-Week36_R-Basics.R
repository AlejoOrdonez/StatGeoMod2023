# Slide 12 - The R way of doing things structure – Functions/Methods + Arguments

sqrt(x = 16)
out = sqrt(x = 16)
out

# Slide 15 - The R way of doing things structure – vector Objects
## Assigning values to an object
a = 10
a
a <- 10
a
10 -> a
a
## Assigning values to an object via operations
a = a + 1
a
b = a * a
b
x = sqrt(b)
x
## Assigning values to an object via concatenation
a = c(4,2,5,10)
a
a = 1:4
a
a = seq(1,10) # How many arguments are passed to this function?
a

# Slide 16 - The R way of doing things structure – matrix objects
A = matrix(data = 0, nrow = 6, ncol = 5) # How many arguments are passed to this function?
A # How many dimensions does this object have

# Slide 17 - The R way of doing things structure – logical Objects
3 < 5
3 > 5
x = 5
x == 5 ## how is this different than the command in line 37?
x != 5 
x = 1:10x < 5

# Slide 20 - The R way of doing things structure – Data frames
## Building a data frame
Data <- data.frame(emp_id = c (1:5),
				  emp_name = c("Rick","Dan",
							   "Michelle","Ryan",
							   "Gary"),
				  salary = c(623.3, 515.2,
							 611.0, 729.0,
							 843.25),
				  start_date = c("2012-01-01",
								 "2013-09-23",
								 "2014-11-15",
								 "2014-05-11",
								 "2015-03-27"))
# Slide 21 - The R way of doing things structure – Data frames
## Print a Data frame
Data ## Notice that the length of each variable/column is the same



# Slide 23 - The R way of doing things structure – Lists
### Build a list
list_data <- list(Text1 ="Red",
				  Text2 ="Green",
				  Numb1 = c(21,32,11),
				  Log1 = TRUE,
				  Nub2 = 51.23)
list_data # Notice that the length of each variable is not the same

# Slide 30 - The R way of doing things Working with data – FINDING your data
setwd(“C:/Project Directory/”)
dir()

# Slide 31 - The R way of doing things Working with data – FINDING your data
setwd(“C:\\Project Directory\\”)
dir()


# Slide 33 - The R way of doing things Working with data – FINDING your data
# write.csv(myData,”updated data.csv”) ### DO NOT RUN!!!
dir()

# Slide 37 - The R way of doing things Packages
a <- available.packages()
head(rownames(a), 3) ## Show the names of the first few packages

# Slide 39 - The R way of doing things loading Packages
library(ggplot2) # might not work if you don't have this one installed
search()[1:6]

	 