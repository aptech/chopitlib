/*
** This file builds and saves the data file
** used to estimate the HOPIT and OP models
*/
new;
cls;
screen on;

// Load raw data and name variables
fname = __FILE_DIR $+ "share_gauss_inc.dat";

// Load full file
xx = loadd(fname);
N = rows(xx);

// Load y variable
y = loadd(fname, "pain");

// Load X variables
x = loadd(fname, "any_cond + grip34_9 + educ_ps");

// Create male variable
male = 1 - loadd(fname, "female");

// Create age range variables
age = loadd(fname, "age");
age6675 = (age .ge 66) .* (age  .le 75);
age76   = (age  .gt 75);

// Concatenate to create full xvars
x = ones(N, 1)~male~x~age6675~age76;
k_all = cols(X);

/*
** Define vignettes by setting V
** This can include v1pain,
** v2pain, and/or v3pain
*/
V = loadd(fname, "v1pain");
kv = cols(v);

// Remove missing values
dataset = packr(y~X~V);

// Extract individual variables
y = dataset[., 1];
x_all = dataset[., 2:cols(dataset)-kv];
V = dataset[., cols(dataset)-kv+1:cols(dataset)];
N = rows(y);

// Variable names
X_vars  = "Constant"$|"Male"$|"Any_cond"$|"Grip34_9"$|"Educ_ps"$|"Age6675"$|"Age76";
V_vars = "V1pain";
Y_vars = "Pain";
all_vars = Y_vars$|X_vars$|V_vars;
if not saved(dataset,"Per_V1_DAT",all_vars);
    errorlog "Write error";stop;
    end;
endif;
