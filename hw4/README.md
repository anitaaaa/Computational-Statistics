This includes the project description solution files. Implemented in R.

Adjustation for hw4:
Dear Class,

Homework 4 is attached.  To complete this homework, you will need to use the Newton-Raphson algorithm for optimization we discussed in class yesterday.  One minor change: You do not need to follow the simulation procedure described in Problem 2 (c) of this assignment if you have trouble with it.  You may use any data generating scheme you wish, provided that the Poisson regression model is correctly specified.  That is, the outcomes need to be Poisson distributed, but you do not need to use these settings in runif or rnorm:

for(i in 1:n) z[i,] = c(1, runif(d-1, -3,1))
beta = rnorm(d,mean=0,sd=1)
