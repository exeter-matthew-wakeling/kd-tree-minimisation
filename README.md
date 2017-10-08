# kd-tree-minimisation
Software to calculate the value of x that produces a given f(x), using a training set.

This software was written by Matthew Wakeling.

## Compiling
Compile the software with the following command:
```
javac PredictFromStats.java
```

## Principles of operation
This software is used to calculate the value of x, when f(x) is known.
Normally, this operation would be performed by analysing f, or by using a minimisation algorithm.
A cost function for minimising is:
```
J(x) = (y - f(x))^T R^(-1) (y - f(x))
```
where y is the known value of f(x), and R is the covariance matrix of the errors in the y vector.

However, some functions are too non-linear for minimisation to work correctly.
This software tackles those difficult functions, by building a database of (R^(-0.5) f(x)) -> x mappings, and then using an efficient multi-dimensional index to find the closest values of f(x) in the database.
The values of x in these closest values is then averaged to produce the resulting estimate of x.

## Usage
First, build a large set of R^(-0.5) f(x) -> x mappings using random values of x.
If some values of x are more likely than others, then it helps to favour these values in the large set.
These mappings should be stored in a tab-separated file, with the x vector and the R^{-0.5) f(x) vector as separate multiple columns.

Secondly, build a similar file containing all the R^(-0.5) f(x) values that you want a value of x for, with the columns in the same positions as the above file. For the columns that hold the x value in the above file, any value can be placed - the number is not used in the prediction, but copied into the output. If the true x values are known (for instance, you are evaluating the performance of the software), they can be put here, and the software will calculate the error in the prediction.

Run the software, as follows:
```
java PredictFromStats -v<number of x columns> <column number> <column number> -c<number of R^(-0.5) f(x) columns> <column number> <column number> -p <number of points to average> --lookup <values to look up> < <training set> > <output>
```
For example:
```
java PredictFromStats -v2 1 2 -c3 3 4 5 -p 200 --lookup unknown.csv <training.csv >output.csv
```
The output file is a tab-separated list, with three columns for each value in x. Of each set of three columns, the first is the value copied from the input file, the second is the calculated value, and the third is the difference between them.

Note that the software detects the number of CPU cores available on the system, and performs that many lookups in parallel. The software requires sufficient RAM to hold the training set, and this can be set using the java option ```-Xmx``` - for example ```java -Xmx3g PredictFromStats...``` for 3GB.

### Other options
+ --limit <number> This forces the software to only use the first <number> rows of the training set, in order to help evaluate how big a training set is required.
+ --retrieved This tells the software to output the coordinates of the training points retrieved from the training set, instead of the prediction made from them. Your list of values to look up should contain just one entry in this case.
