{smcl}
{*  16 July 2015}{...}
{cmd:help Rpovmap}
{hline}


{title:Title}

    Wrapper for running the R package ‘povmap’ for small area estimation from within Stata

	
{title:Syntax}

{p 8 17 2}
{opt Rpovmap} {depvar} {indepvars} {cmd:,}
{cmdab:smp_data(}{it:string}{cmd:)}
{cmdab:pop_data(}{it:string}{cmd:)}
{cmdab:smp_domains(}{it:string}{cmd:)}
{cmdab:pop_domains(}{it:string}{cmd:)}
[{it:options}]


{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt pop_data(string)}} specifies the file name of the population data file. Please include the .dta extension. If no path is given, the current directory is assumed. {p_end} 
{synopt:{opt pop_domains(string)}} specifies the name of the variable that indicates domains in the population data. {p_end} 
{synopt:{opt smp_data(string)}} specifies the file name of the sample data file. If no path is given, the current directory is assumed. {p_end} 
{synopt:{opt smp_domains(string)}} specifies the name of the variable that indicates domains in the sample data. {p_end} 

{title:Options}
{synopt:{opt l(integer)}} a number determining the number of Monte-Carlo simulations that must be at least 1. Defaults to 50 {p_end}
{synopt:{opt threshold(string)}} a number defining the poverty line threshold for calculating headcount rates. Defaults to 60% of the median of the dependent variable {p_end}
{synopt:{opt transformation(string)}} specifies the transformation to be applied to the dependent variable. Options are "no","log","box.cox","dual","log.shift",and "ordernorm". Defaults to "box.cox"  {p_end} 
{synopt:{opt interval(string)}} a string equal to 'default' or a numeric vector containing a lower and upper limit determining an interval for the estimation of the optimal parameter. Defaults to "c(-1,2)"  {p_end} 
{synopt:{opt MSE(string)}} if TRUE, MSE estimates using a parametric bootstrap approach are calculated. Defaults to FALSE {p_end} 
{synopt:{opt b(integer)}} A number determining the number of bootstrap populations in the parametric bootstrap approach used in the MSE estimation. The number must be greater than 1. Defaults to 50. For practical applications, values larger than 200 are recommended {p_end} 
{synopt:{opt seed(integer)}} an integer to set the seed for the random number generator. For the usage of random number generation. Defaults to 123 {p_end} 
{synopt:{opt boot_type(string)}} character string to choose between different MSE estimation procedures. The two options are "parametric" and a semi-parametric "wild" bootstrap. Defaults to "parametric". {p_end} 
{synopt:{opt cpus(integer)}} an integer determining the number of cores used when estimating the MSE. Defaults to 1. {p_end} 
{synopt:{opt na_rm(string)}} if TRUE, observations with NA values are deleted from the population and sample data. For the EBP procedure complete observations are required. Defaults to FALSE. {p_end} 
{synopt:{opt weights(string)}} specifies the name of the weight variable contained in the sample data. Default is no weights. {p_end} 
{synopt:{opt pop_weights(string)}} specifies the name of the weight variable contained in the population data to use when aggregating estimates. Default is no weights.  {p_end} 
{synopt:{opt aggregate_to(string)}} a character string containing the name of a variable from the population data to calculate estimates for. Default is the variable specified in the smp_domain option.  {p_end} 
{synopt:{opt weights_type(string)}} a character string. Three different options for survey weights are available (i) EBP under informative sampling from Guadarrama et al. (2018) ("Guadarrama"); (ii) considering survey weights by using the weighting options of nlme from Pinheiro and Bates (2023) when estimating the linear mixed model ("nlme"); (iii) considering survey weights by using the weighting options of nlme and also using the specified weights when determining the optimal transformation parameter lambda ("nlme_lambda"). Defaults to "Guadarrama". {p_end} 
{synopt:{benchark(string)}} An R vector containing the names of the indicators to be benchmarked, i,e. a value of `"c("Mean","Headcount")"' to benchmark both the mean and the headcounta. Unlike the R package, only internal benchmarking using the sample data is available. Defaults to `"c("Head_Count")"' {p_end}
{synopt:{benchark_type(string)}}  A character indicating the type of benchmarking. The two options are (i) Raking ("raking") and (ii) Ratio adjustment ("ratio"). Defaults to "ratio"' {p_end}
{synopt:{benchark_level(string)}}  A character specifying the variable in the sample and population data indicating the level at which the benchmarking is performed. This option must be specified if the benchmark option is selected. {p_end}
{synopt:{benchark_weights(string)}}  A character specifying the variable in the sample used to weight the sample data when calculating the values to bencmhark to. Defaults to the variable specified in the "weights" option {p_end}
{synopt:{nlme_maxiter(integer)}} an integer indicating the maximum number of iterations the lme function from package nlme will carry out before exiting with a non-convergence error message. Defaults to 1000. 
{synopt:{nlme_tolerance(real)}}  an integer indicating the tolerance criterium for the the lme function from package nlme. Defaults to 1e-6
{synopt:{rescale_weights(string)}}  A logical indicating if the weights in the sample data are scaled. If FALSE (default), the sample weights do not change. When TRUE, the sample weights are rescaled such that the average weight is 1 within each target domain. {p_end}




{title:Description}

{pstd}
Allows small area estimates using Empirical Best Predictor models using the R package ‘povmap’ from within Stata. 


{title:Dependencies}

{pstd}
{hi:R} -- The R software system must be installed on the user's system in order for {hi:Rpovmap} to run. You can download R from https://www.r-project.org)
{hi:R} -- Rpovmap requires the haven and povmap packages to be installed in R. Executing install.packages("haven") and install.packages("povmap") in R will install these packages from CRAN 

{pstd}
{hi: rscript} -- Additionally, the Stata package rscript is used to run an R source file from within Stata. You can type "ssc install rscript" to install the rscript package

{title:Examples}

{pstd}
{hi:Example}: Estimating a model 

{pstd}
{stata sysuse auto, replace} <---- First, load everyone's favorite dataset.

{pstd}
Use {hi:Rlasso} to estimate a Lasso model using 5 fold cross validation. {p_end}

{phang}
{stata gen lnprice = log(price)}

{phang}
{stata Rlasso lnprice mpg rep78 headroom trunk weight length turn displacement gear_ratio foreign, nfolds(5)}

{pstd}
{hi:Example 2}: Predicting low birth weight using a logit model with ridge penalty, reporting coefficients at the MSE minimizing lambda value.

{phang}
{stata "use http://www.stata-press.com/data/r13/lbw.dta, clear"}

{phang}
{stata Rlasso low age lwt race smoke ptl ht ui ftv, model(logit) lambda(min) penalty(0)} 


{title:Saved Results}

{pstd} 
{hi:Rlasso} displays all results to the results window. 


{title:Installation Instructions}

{pstd}
{hi:1.} Ensure you have R executable installed. If you don't have this installed, you can download and install the software from https://www.r-project.org/. 

{pstd}
{hi:2.} Install the Stata dependency Rsource in Stata. {stata "ssc install Rsource"}

{pstd}
{hi:3.} If you have a windows machine, the location of the Rterm executable file may need to be manually set using the global {hi:Rterm_path} or the option {hi: rpath}. First find the location of the Rterm executable file on your local machine. On our test PC this is installed in the folder "C:\Program Files\R\R-3.1.2\bin\i386\Rterm.exe". Set the path to this executable using Rterm_path using "global Rterm_path [path]" or via the {hi: rpath} option. 


    
	
{title:Authors}

{pstd} 
Jonathan Hersh, Poverty Global Practice, World Bank, Washington, D.C. 
Ana Areias, Poverty Global Practice, World Bank, Washington, D.C. 

{pstd} 
Email {browse   "mailto:jhersh@worldbank.org":jhersh@worldbank.org} 
Email {browse   "mailto:ana.areias@gmail.com":ana.areias@gmail.com} 


{title:References}

{pstd} 
Friedman, J., Hastie, T. and Tibshirani, R. (2008) {it:Regularization Paths for Generalized Linear Models via Coordinate Descent}

{pstd} 
Rsource filelist

