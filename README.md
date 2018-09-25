# extremes

This repo contains the code needed to calculate and plot the results of  "Estimating Changes in Temperature Distributions in a Large Ensemble of Climate Simulations Using Quantile Regression" (10.1175/JCLI-D-17-0782.1) Journal of Climate, vol. 31, no. 20. 

There are a number of key files:
1. boot2d.R - This is the file where the main computation is performed, i.e. the quantiles are calculated for North America.
boot2d-mpi.R Is the same thing but using a large cluster package `pbdMPI`. The ...-mx/-mn are for daily maximum and daily minimum temperatures.
Note that all these files contain hardcoded file paths which is a result of the proprietary nature of the raw data. 
That is why these files will not run on another local machine.
These files should rather be used as guidelines on how to calculate the quantiles of interest. 
The main line which executes the quantile fit is the following
```coef[,i] = rq.fit.pfn(X, y=ts_data, tau=q[i], max.bad.fixup=20)$coefficients```, which relies on the `quantreg` package.
The `helper.R` file contains a lot of helper functions used in plotting etc.

2. mainBig-mpi.R performs the same calculation as boot2d.R but considers a larger area.

3. CV.R Is the model selection module where I calculate which model is the best one, where the main calculating 
loop is on line 104 where we cross-validate one pixel. The rest of the module does some plotting.

4. The paper.R file contains the code for plotting most of the figures in the paper.