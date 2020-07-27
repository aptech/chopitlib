new;
cls;
screen on;
clearg jj, kx_struct, kx_mu, kv, kx_all, capJ, xbar, X_bar_struct, X_bar_mu, p_0, exp_mu, check_con, cond_Y_N, sig_restrict;


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

		/* PROGRAM TO ESTIMATE VIGNETTES HOPIT MODEL */

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


library cmlmt;

/* SET DIRECTORIES */

#include "Vignettes_LM_test.prc";

/* LOAD RAW DATA AND NAME VARIABELS */

infile = "share_gauss_inc";

open fh = ^infile varindxi;
vnames2  = getname(infile);
	   
xx 	= readr(fh,rowsf(fh));
N	= rows(xx);

/* DEFINE DEPENDENT VARIABLE AND ID VARIABLE */

y			= xx[., ipain];
let yvars   = PAIN;

idvec   	= seqa(1, 1, n);
male		= 1-xx[., ifemale];
	
/* DEFINE COMMON VARIABLES: THESE WILL FEATURE IN ALL PARTS OF THE MODEL */ //Note: these must be same as vignettes vars above

/* *** user defined set of explanatory variables *** */

age6675 = (xx[.,iage] .ge 66) .* (xx[.,iage] .le 75);
age76   = (xx[.,iage] .gt 75);

X 		    = ones(N, 1)~male~xx[., iany_cond igrip34_9 ieduc_ps]~age6675~age76;
let X_vars  = constant male any_cond grip34_9 educ_ps age6675 age76;
k_all		= cols(X);

	
/* DEFINE VIGNETTES */

V1		= xx[., iv1pain];
V2		= xx[., iv2pain];
V3		= xx[., iv3pain];

/* *** stack what vignettes to be used together in "V" and name accordingly *** */
    
//V		= V1~V2~V3; 
//V		= V1~V3; 
//V		= V1~V2; 
//V		= V2~V3; 
//V		= V1;
//V		= V2;
V		= V1;

//let V_vars  = vig1 vig2 vig3;
//let V_vars  = vig2 vig3;
//let V_vars  = vig1 vig2;
let V_vars  = v1pain;

kv			= cols(V);
		
/* DELETE ANY MISSING VALUES */

ydata      	= y~X~V;
ydata       = packr(ydata);
y 	        = ydata[.,1];
X_all      	= ydata[.,2:cols(ydata)-kv];
av_X      	= meanc(X)';
V			= ydata[.,cols(ydata)-kv+1:cols(ydata)];
N       	= rows(y);

/* CONSTRUCT INDICATOR VARIABLE FOR j=J */

jrep = minc(y);
    
d_ij	= zeros(rows(y),1);
V1_ij   = zeros(rows(y),1);
   
do while jrep lt maxc(y)+1;
        
	d_ij = d_ij~(y .eq jrep);
    V1_ij = V1_ij~(V .eq jrep);
        
jrep=jrep+1;
endo;

d_ij    = d_ij[.,2:cols(d_ij)];
V1_ij   = V1_ij[.,2:cols(V1_ij)];    

phat    = meanc(d_ij);
phat_V1 = meanc(V1_ij);

ydata 	= y~X_all~V; 

all_vars    	= yvars|X_vars|V_vars;
    
if not saved(ydata,"Per_V1_DAT",all_vars);
    errorlog "Write error";stop;
    end;
endif;
    
{ vnam, mean, var, std, min, max, valid, mis } = dstat("Per_V1_DAT", all_vars);
    

/* SET EXP_MU=1 FOR EXP(.) THRESHOLDS; =0 FOR LINEAR ONES; TEST ONLY COMPUTED FOR EXP_MU = 1 */

exp_mu = 1;


/* check_con = 1 re-estimate the HOPIT model with solved start values to ensure global maximum. If HOPIT model doesn't look to have converged properly, turn this on */

check_con = 0;

/* condition LM tests? If the test fails to compute, try conditioning it with a very small number here */

cond_Y_N = 0;


/* PRINT ITERATIONS? NO = 0; YES = 1*/

print_it = 1;

/* SET TOLERANCE LEVEL (-1 FOR DEFAULTS) */

tol = -1;

/* RESTRICT ALL SIGMAS=1? (see paper) */

sig_restrict = 1;


/* ENTER OPTION START VALUES: to be taken from the Gauss CMLMT output screen, in the exact same order */

let start_OP = 0
;

let start_HOP =   0;

/* REQUEST APPROPRIATE MODEL: OP = just Ordered Probit; HOPIT estimates the CHOPIT vignettes model */

let model 	= OP;
let model   = HOPIT;

let mtype = OP HOPIT;
mtype = indcv(model,mtype); 


/* CALL PROCEDURE */
		
	{OLS_B, 
		OP_B, ME_OP, ME_OPse, V_OP, L_OP, IC_OP, pordp, y_star_OP, 
		HOPIT_B, me_HOPIT, me_HOPITse, V_HOPIT, L_HOPIT, IC_HOPIT, p_HOPIT, y_star_HOPIT,   
		IC_ALL, IC_HOPITS, LM_test}

		= CHOPITS_LM(0,"Per_V1_DAT", start_OP, start_HOP, print_it, tol, model);	

reg_vars	= X_vars;
OP_vars		= X_vars[2:rows(reg_vars)];

/* PRINT RESULTS */

if exp_mu eq 1;
    output file = "Perachi_J4_ammended_exp_V3.out" reset;
else;
     output file = "Perachi_J4_linear_exp.out" reset;
endif;

output on;

format 15,4;

?;
"               SELF-ASSESSED HEALTH (C)CHOPIT MODEL";
?;
"DATE IS";;datestring(0);
?;
if exp_mu eq 0;
	"LINEAR THRESHOLDS FOR HOPIT VARIANTS";
elseif exp_mu eq 1;
	"EXPONENTIAL THRESHOLDS FOR HOPIT VARIANTS";
endif;	

"DESCRIPTIVE STATISTICS: SAH DATA";?;

{ vnam, mean, var, std, min, max, valid, mis } = dstat("Per_V1_DAT",all_vars);
?;

" No. of vignettes used";;kv;?;
"Vignette names";;$V_vars[1:kv];
?;

"Sample Proportions of Y";?;
phat;
?;
	
if exp_mu eq 1;

	"LM TEST FOR RC AND VE";
	"Test, df's and p-value (reject null HOPIT model for large test values and/or small p-value)";?;
	LM_test[1,.];
	?;
	"LM TEST FOR VE";
	"Test, df's and p-value (reject null HOPIT model for large test values and/or small p-value)";?;
	LM_test[2,.];
	?;
	"LM TEST FOR RC";
	"Test, df's and p-value (reject null HOPIT model for large test values and/or small p-value)";?;
	LM_test[3,.];
	
endif;	

?;
"AVERAGE PROBABILITIES";
?;

if mtype eq 1;
	meanc(pordp);		
elseif mtype eq 2;		
	meanc(pordp)~meanc(p_HOPIT);
endif;

?;
"SUMMARIES OF INFORMATION CRITERIA (BIC, AIC, CAIC, HQIC)";?;
format 15,8;
"OP REGRESION";miss(IC_OP,1000000);?;
if mtype ge 2;
	"HOPIT";miss(IC_HOPIT,1000000);?;
endif;

?;
if mtype ge 2;
	"No. of vignettes used";;kv;?;
	"Vignette names";;$V_vars;
endif;

print;
print;
"		SIMPLE REGRESSION";
format 15,4;

print;
"ESTIMATED PARAMETERS: COEFFICIENT S.E. T-STAT";
print;

for i (1,kx_struct,1);
	$reg_vars[i];;OLS_B[i,.];
endfor;
print;

let mu 		= mu;
mu			= mu .* ones(capJ-1,1);

OP_vars 	= reg_vars[2:rows(reg_vars)]|mu;
$OP_vars;

print;
print;
"		SIMPLE ORDERED PROBIT REGRESSION";
format 15,6;

print;
"ESTIMATED PARAMETERS: COEFFICIENT S.E. T-STAT";
print;

for i (1,kx_struct-1+capJ-1,1);
	$OP_vars[i];;OP_B[i,.];
endfor;
print;

print;
"ESTIMATED PARTIAL EFFECTS (j=1,...,J; standard errors underneath)";
print;

for i (1,kx_struct-1,1);
	$OP_vars[i];;ME_OP[i,.];
	"		";;ME_OPse[i,.];
endfor;
print;

if mtype ge 2;

	gama		= HOPIT_B[kx_struct+kv+1:kx_struct+kv+(kx_mu*(capJ-1)),1];
	gama	    = reshape(gama,capJ-1,kx_mu)';  
	gama		= vec(gama);
	
	gama_se		= HOPIT_B[kx_struct+kv+1:kx_struct+kv+(kx_mu*(capJ-1)),2];
	gama_se		= reshape(gama_se,capJ-1,kx_mu)';  
	gama_se		= vec(gama_se);
	
	gama_t		= HOPIT_B[kx_struct+kv+1:kx_struct+kv+(kx_mu*(capJ-1)),3];
	gama_t		= reshape(gama_t,capJ-1,kx_mu)';  
	gama_t		= vec(gama_t);
	
	gama_all	= gama~gama_se~gama_t;

	print;
	print;
	"		HOPIT VIGNETTES MODEL";
	format 15,6;

    print;
    "LogL";;L_HOPIT;
    print;

	print;
	"ESTIMATED PARAMETERS: COEFFICIENT S.E. T-STAT";
	print;
	"Structural Parameters"; 
	print;

	for i (1,kx_struct,1);
		$reg_vars[i];;HOPIT_B[i,.];
	endfor;
	print;

	"Vignettes constants";
	print;
	
    for i (1,kv,1);
            $V_vars[i];;HOPIT_B[kx_struct+i,.];
    endfor;
    print;
        
	"Boundary equations (j=0,...,J)";
	print;
	
	for j (1,capJ-1,1);
	
		for i (1,kx_mu,1);
			$reg_vars[i];;gama_all[((j-1)*kx_mu)+i,.];
		endfor;
		print;
		
	endfor;

    if rows(HOPIT_B) gt kx_struct+kv+((capJ-1)*kx_mu);
        
        "	 1/V(sd)";;HOPIT_B[kx_struct+kv+((capJ-1)*kx_mu)+1,.];
        print;
        
    endif;
    
	print;
	"ESTIMATED PARTIAL EFFECTS (j=1,...,J; standard errors underneath)";
	print;

	ME_HOPIT	= delif(ME_HOPIT,ME_HOPIT[.,1] .eq 1); 
	ME_HOPITse 	= delif(ME_HOPITse,ME_HOPITse[.,1] .eq 1); 
	
	for i (1,kx_struct,1);
		$X_vars[i];;ME_HOPIT[i,.];
		"		";;ME_HOPITse[i,.];
	endfor;
	print;

endif;


output off;
closeall;

stop;

