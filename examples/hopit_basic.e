new;
library cmlmt, chopitlib;

// Declare filename for running code
fname = __FILE_DIR $+ "Per_V1_DAT.dat";

// Run descriptive statistics on dataset
dstatmt(fname);

// Set up variable names 
// for data groups
v_formula = "v1pain";
x_formula = "Constant + Male + any_cond + grip34_9 + educ_ps + age6675 + age76";
y_formula = "pain";

// Call Procedure
struct chopitOut cOut;
cOut = CHOPITS_LM(fname, v_formula, y_formula, x_formula);	
