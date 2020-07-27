new;
library cmlmt, chopitlib;

// Declare filename for running code
fname = __FILE_DIR $+ "Per_V1_DAT.dat";

// Run descriptive statistics on dataset
dstatmt(fname);

/*
** This section sets up the parameters for
** running the CHOPIT model
*/
struct chopitControl cCtl;

// Set up thresholds using exp_mu
//     0       for linear ones   
//     1       exponential thresholds
//             (test only computed ones)
cCtl.exp_mu = 1;

// Turn on check for convergence
// Use this if you have concerns about
// convergence
//    0        no check
//    1        re-estimate HOPIT model with solved start
//             values to check for global convergence
cCtl.check_con = 0;

// Condition LM tests
// If the test fails to compute, try conditioning 
// it with a very small number here 
cCtl.cond_Y_N = 0;

// Printing of iterations 
//    0        No 
//    1        Yes 
cCtl.print_it = 1;

// Set tolerance level
// -1 for defaults
cCtl.tol = -1;

// Restrict all sigmas 
cCtl.sig_restrict = 1;

// Enter option start values :
// To be taken from the Gauss CMLMT 
// output screen, in the exact same order
cCtl.start_OP = 0;
cCtl.start_HOP = 0;

// Request appropriate model:
// OP     = just ordered probit
// HOPIT  = estimates the CHOPIT vignettes models
cCtl.model = "HOPIT";

// Set up variable names 
// for data groups
cCtl.v_varnames = "v1pain";
cCtl.x_varnames = "Constant"$|"Male"$|"any_cond"$|"grip34_9"$|"educ_ps"$|"age6675"$|"age76";
cCtl.y_varnames = "pain";

// Call Procedure
struct chopitOut cOut;
cOut = CHOPITS_LM(0, fname, cCtl);	
