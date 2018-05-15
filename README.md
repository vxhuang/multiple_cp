# multiple_cp
Multiple changepoint detection

example.R shows two examples. The first is a times series with 2 changepoints, and the second is a times series with 1 changepoints. This algorithm can differentiate them pretty accurately. 

amtc.R is an algorithm designed to detect changepoints from time series with at most two changepoints. 

Suppose that two changepoints are known to occur at time $$c_1$$,$c_2$, we can use genetic algorithm to discover the optimal solution.
Mutation: add or subtract c_1 or c_2 by a random small amount;
Crossover: exchange the c_1 and c_2 from two selected solutions.

In my trial, this algorithm converges really fast (typically in ~5 iterations) but may fall in sub-optimal solution. To achieve a better performance, we need a relatively high volume of initial choices. 

We should always notice that the choices of the changepoints should not be too close to the boundaries or too close to each other. 
