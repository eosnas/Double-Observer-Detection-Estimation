# Double-Observer-Detection-Estimation
Functions and code for estimating detection rate from two-occasion observation data (double observer).  

## Functions for making matches
Some functions for determining a match in double observer data from front seat and back seat aerial observations used for estimating detection probability. These function would on files structure that might be different than currently used (i.e., the quality control process and column names of the input data may have changed).

You will need to read the code to understand the function, there is only minimal and insufficient documentation and comments within in the function code. A brief description of each file is :

**DubMatch.R** is a minial function that makes the most strict assumption about matches. There is no identification error and unit = 'open' and 'single' are treated as matches. The time window for matches is set as input to the function, default is appropriate for a plane moving at about 45 meters per second. This function should be tested to ensure it works as intended.

**DubMatch2.R** is a liberal matching scheme where different units can be use by front or back seat observes and still a match is allow. See function for all combinations but include units with multiple 'singles' or 'opens' seen in the time window, and also allow 'singles' and 'floked drakes' to be matched for various group sizes. Allows 'singles' to be matched to 'pairs' and small flocks within 3 to 5 to match even if number is not exact. Any group > 5 is a match if the unit combination is reasonable. See code.

**findNear.R** is a function that find rows of a data frame that are close in time or space (i.e., observation within a slinding window).

**formatDDdata.R** is a function that cleans up double observer data for use in the other function. Structure of this old data is different than current observation data, so this function might not be relavent for double observer data collected after 2016. However, this function serve to document some of the quality control isses in the 2015 and 2016 data.

**matching_matrix.csv** is a file that defines a matrix of front seat and back seat combinations of unit and group size can reasonably be considered matches. It might be better to program the match function so that a matrix like this could be input to define matches, rather than a large series of if then statments. This would make the function more general and probably a lot easier to program.
