library cmlmt;
#include cmlmt.sdf;

/* spline_RP_OP_MT.PRC: RANDOM PARAMETERS OP SPLINE MODEL */

cls;
trap 0;

proc 20 = CHOPITS_LM(dset_V, dset, start_OP, start_HOP, print_it, tol, model);
    local infile_V, infile, fh_V, fh, varsz, y, yVX, V, X, X_struct, X_mu, X_inf, n, zcheck, zchkvec, cdo, id, mtype, V_x, VX;
    local y_star_OP, y_star_HOPIT, y_star_MIOP, y_star_TOP, y_star_MIOP_HOP, y_star_TOP_HOP;
    local OLS_B, OP_B, me_OP, me_OPse, V_OP, L_OP, IC_OP,
        HOPIT_B, me_HOPIT, me_HOPITse, V_HOPIT, L_HOPIT, IC_HOPIT,
        MIOP_B, me_MIOP, me_MIOPse, V_MIOP, L_MIOP, IC_MIOP;
    local	V1, V2, V3, V4, V5, d_ij, phat, v_ij1, v_ij2, v_ij3, v_ij4, v_ij5, nam, m, stb, vc_ols, std, sig, cx, rsq, resid, dbw;
    local startmu1, startmu, startrep, mu_temp, beta_Y, c, op_mu, retOPb, t_OP, et_OP, beta_struct, OPbse, zbetahat, mu_mat, mu_xbeta, p_ij;
    local pordp, check, me_model, pOPxbar, pOP_se, V_cons, sig_V_inv, gama, K_HOPIT, retHOPIT, et_HOPIT, HOPIT_se, t_HOPIT;
    local zgama, mu_ij, xbeta, p_HOPIT, p_HOPIT_xbar, p_HOPITse, y_inf, sel, beta_inf, beta_inf_struct, retMIOP, et_MIOP;
    local mu_MIOP, MIOP_se, t_MIOP, xb_inf, xb_inf_struct, xb_struct, p_MIOP, p_MIOP_prior, p_no_inf, p_no_inf_prior, V_GOP_V2;
    local p_aug0, p_aug1, p0_MIOP, p_inf_MIOP, pLC2_0, pLC2_1, p_inf_prior, L0, L1, L, LC_predict, pPOST_MIOP, p_MIOP_xbar, p_MIOPse;
    local m_i, Vuong_OP_MIOP, lnLiHOPIT, GOP_1_se;
    local TOP_B, me_TOP, me_TOPse, V_TOP, L_TOP, IC_TOP, beta_struct_TOP, mu_TOP, tempered1_b, tempered2_b;
    local tempered1_mu, tempered2_mu, tempered3_b, tempered3_mu, retTOP, et_TOP, mu1, mu2, mu3, TOP_se, t_TOP;
    local p_TOP, p_TOP_prior, xb_temp1, p_temp1, xb_temp2, p_temp2, mu_mat1, mu1_xbeta, mu_mat2, mu2_xbeta, mu_mat3, mu3_xbeta, xb_temp3, p_temp3;
    local p_TOP_xbar, p_TOPse, Vuong_OP_TOP, Vuong_MIOP_TOP, IC_ALL, IC_HOPITS;
    local k_MIOP_HOP, L_MIOP_HOP, V_MIOP_HOP, retMIOP_HOP, et_MIOP_HOP, MIOP_HOP_B, MIOP_HOP_se, t_MIOP_HOP, p_MIOP_HOP, p_MIOP_HOP_prior;
    local p0_MIOP_HOP, p_inf_MIOP_HOP, pPOST_MIOP_HOP, me_MIOP_HOP, me_MIOP_HOPse, p_MIOP_HOP_xbar, p_MIOP_HOPse, Vuong_HOP_MIOP_HOP, IC_MIOP_HOP;
    local TOP_HOP_B, me_TOP_HOP, me_TOP_HOPse, V_TOP_HOP, L_TOP_HOP, IC_TOP_HOP, p_TOP_HOP, p_TOP_HOP_prior, Vuong_HOP_TOP_HOP;
    local Vuong_MIOP_HOP_TOP_HOP, retTOP_HOP, et_TOP_HOP, TOP_HOP_se, t_TOP_HOP, p_TOP_HOP_xbar, p_TOP_HOPse;
    local gama_temp, gama0, alpha0_free, beta0_free, delta0_free, k_GOP, L_GOP, V_sig;
    local gama1_free, delta1_free, gama2_free, delta2_free, gama3_free, delta3_free, gama4_free, delta4_free, theta0_free, theta1_free, theta2_free, V_sigma;
    local psy_start, R_gama;
    local beta_V1, gama0_V1, gama_V1, beta_V2, gama0_V2, gama_V2, beta_V3, gama0_V3, gama_V3, k_V_GOP, L_V_GOP;
    local gama0_HOP, gama_HOP, beta_HOP, LR_test, df, LR_p;
    local gama0_MIOP, gama_MIOP, gama0_TOP, gama_TOP, k_TOP_HOP, se_gama;
    local L_GOP_V1, L_GOP_V2, L_GOP_V3, L_un, no_param_un, k_V1, k_V2, k_V3, df_beta, df_gama0, df_gama, V_GOP, GOP_se, params;
    local lnL_SAH, lnL_V1, beta_LR, gama_LR, mu_ijV, mu_xbetaV1, p_ijV1, lnLiV1;
    local R, R_b_primary, R_gamma_plus, R_gamma_minus, q, free_params, Wald, params_V1, R_b_vignettes, V_all, V_GOP_V1;
    local V_q, Wald_p, Wald_test, k_GOP_JOINT, V_GOP_JOINT, GOP_JOINT_se, L_GOP_JOINT, df_beta0, df_V_beta, R_beta0;
    local count_hop, count_gop, score, LM_test, grad_beta, grad_g1, grad_g2, grad_g3, grad_g4;
    local xbeta_V1, mu_xbeta_V1, p_ij_V1, grad_beta_V1, grad_g1_V1, grad_g2_V1, grad_g3_V1, grad_g4_V1, Score_anal, bhhh_anal;
    local params_V2, grad_beta_V2, grad_g1_V2, grad_g2_V2, grad_g3_V2, grad_g4_V2, xbeta_V2, mu_xbeta_V2, p_ij_V2, R_1V, R_2V, bhhh_max;
    local beta_check, xbeta_V3, mu_xbeta_V3, p_ij_V3, grad_beta_V3, grad_g1_V3, grad_g2_V3, grad_g3_V3, grad_g4_V3;
    local Score_anal_RC, Score_anal_VE, df_RC, bhhh_anal_RC, df_VE, bhhh_anal_VE, V_RC, V_VE, LM_test_RC, LM_test_VE;
    local mu_xbetaV2, p_ijV2, mu_xbetaV3, p_ijV3, fail_linear, chol_bhhh, chol_bhhh_RC, chol_bhhh_VE, V_inv, V_inv_RC, V_inv_VE;
    local condition, cond_bhhh, cond_bhhh_VE, LM_cond, LM_cond_VE, HOPIT_grad, cond_bhhh_RC, LM_cond_RC;
    local av_boundaries, V_gamma, G, p0, var_boundaries, se_boundaries, V_mu0, G0, V_gamma0;
    local mu_ij_exp_new, mu_ij_linear, mu0_star, mu1_star, mu2_star, params_V3, grad_theta_V1, grad_theta_V2, grad_theta_V3, grad_theta;
    local HOPIT_fail, retGOP_JOINT, GOP_JOINT_fail, ret_GOP, ret_GOP_V1, Hessian_joint, Hessian_joint_i, Hessian_VE, Hessian_RC;
    local state, rand, LR_Yes_No;
    
    clearg index0, index1;
    
    clear OLS_B,
        OP_B, me_OP, me_OPse, V_OP, L_OP, IC_OP, pordp, y_star_OP,
        HOPIT_B, me_HOPIT, me_HOPITse, V_HOPIT, L_HOPIT, IC_HOPIT, p_HOPIT, y_star_HOPIT,
        MIOP_B, me_MIOP, me_MIOPse, V_MIOP, L_MIOP, IC_MIOP, p_MIOP, y_star_MIOP, p_MIOP_prior,
        TOP_B, me_TOP, me_TOPse, V_TOP, L_TOP, IC_TOP, p_TOP, y_star_TOP, p_TOP_prior,
        MIOP_HOP_B, me_MIOP_HOP, me_MIOP_HOPse, V_MIOP_HOP, L_MIOP_HOP, IC_MIOP_HOP, p_MIOP_HOP, y_star_MIOP_HOP, p_MIOP_HOP_prior,
        TOP_HOP_B, me_TOP_HOP, me_TOP_HOPse, V_TOP_HOP, L_TOP_HOP, IC_TOP_HOP, p_TOP_HOP, y_star_TOP_HOP, p_TOP_HOP_prior,
        IC_ALL, IC_HOPITS,
        Vuong_OP_MIOP, Vuong_OP_TOP, Vuong_MIOP_TOP, Vuong_HOP_MIOP_HOP, Vuong_HOP_TOP_HOP, Vuong_MIOP_HOP_TOP_HOP, LR_test, Wald_test, LM_test;
    
    clear k_V1, k_V2, k_V3;
    
    LR_Yes_No = 2;      //this was an option before, but turn off
    
    LM_test = zeros(3, 3);
    
    /* WHICH MODELS? */
    
    let mtype = OP HOPIT;
    mtype = indcv(model, mtype);
    
    /* READ IN COMPOSITE DATASET SET */    
    infile 			= dset;
    open fh 		= ^infile varindxi;
    
    yVX   			= readr(fh, rowsf(fh));
    yVX    			= packr(yVX);
    y 	    		= yVX[., 1];
    
    X     			= yVX[., 2:cols(yVX)-kv];
    V				= yVX[., cols(yVX)-kv+1:cols(yVX)];
    
    kx_all			= cols(X);
    X_struct		= X;
    kx_struct       = cols(X_struct);
    X_mu			= X;
    kx_mu           = cols(X_mu);
    
    //cols(X_mu)~cols(X_struct);stop;
    
    X_bar_struct	= meanc(X_struct);
    X_bar_mu		= meanc(X_mu);
    
    y       = y - minc(y);
    capJ    = maxc(y) - minc(y) + 1;
    N       = rows(y);
    fh      = close(fh);
    
    if capJ lt 3;
        "ERROR: too few categories for estimation (3 minimum). # is";;
        capJ;
        stop;
    endif;
    
    V1		= V[.,1] - minc(V[.,1]);
    if kv ge 2;
        V2		= V[.,2] - minc(V[.,2]);
    endif;
    if kv ge 3;
        V3		= V[.,3] - minc(V[.,3]);
    endif;
    if kv ge 4;
        V4		= V[.,4] - minc(V[.,4]);
    endif;
    if kv ge 5;
        V5		= V[.,5] - minc(V[.,5]);
    endif;
    
    /* CONSTRUCT INDICATOR VARIABLE FOR j=J */
    
    d_ij    = zeros(N,1);
    
    for jrep (0,capJ-1,1);
        d_ij = d_ij~(y .eq jrep);
    endfor;
    
    d_ij = d_ij[.,2:capJ+1];
    phat = meanc(d_ij);
    
    v_ij1 = {};
    
    for jrep (minc(V1),maxc(V1),1);
        v_ij1 = v_ij1~(V1 .eq jrep);
    endfor;
    
    if kv ge 2;
        
        v_ij2 = {};
        
        for jrep (minc(V2),maxc(V2),1);
            v_ij2 = v_ij2~(V2 .eq jrep);
        endfor;
        
    endif;
    
    if kv ge 3;
        
        v_ij3 = {};
        
        for jrep (minc(V3),maxc(V3),1);
            v_ij3 = v_ij3~(V3 .eq jrep);
        endfor;
        
    endif;
    
    if kv ge 4;
        
        v_ij4 = {};
        
        for jrep (minc(V4),maxc(V4),1);
            v_ij4 = v_ij4~(V4 .eq jrep);
        endfor;
        
    endif;
    
    if kv ge 5;
        
        v_ij5 = {};
        
        for jrep (minc(V5),maxc(V5),1);
            v_ij5 = v_ij5~(V5 .eq jrep);
        endfor;
        
    endif;
    
    /* ESTIMATE SIMPLE OLS */
    
    OLS_B	= (inv(X_struct'X_struct)*X_struct'y);
    
    
    /* ESTIMATE ORDERED PROBIT MODEL */
    
    if rows(start_OP) eq 1;
        
        startmu1    = meanc(d_ij);
        startmu     = startmu1[1];
        
        for startrep (2,capJ-1,1);
            mu_temp = ((startmu[startrep-1] + startmu1[startrep]));
            startmu  = startmu|mu_temp;
        endfor;
        
        startmu 	= cdfn(startmu);
        
        beta_Y  		= OLS_B[2:rows(OLS_B),1];
        
    else;
        
        beta_Y  	= start_OP[1:kx_struct-1];
        startmu     = start_OP[kx_struct:rows(start_OP)];
        
    endif;
    
    struct PV pORDERED;
    pORDERED = pvCreate;
    pORDERED = pvPack(pORDERED,beta_Y,"beta_Y");
    pORDERED = pvPack(pORDERED,startmu,"mu");
    
    struct DS dORDERED;
    dORDERED = reshape(dsCreate,2,1);
    dORDERED[1].DataMatrix = d_ij;
    dORDERED[2].DataMatrix = X_struct[.,2:cols(X_struct)];
    
    struct cmlmtControl cORDERED;
    cORDERED        = cmlmtControlCreate;
    cORDERED.title  = "Ordered Probit Model";
    
    c           = ((eye(capJ-2) .*-1)~zeros(capJ-2,1)) + (zeros(capJ-2,1)~eye(capJ-2));
    cORDERED.C  = zeros(capJ-2,kx_struct-1)~c;
    cORDERED.D  = zeros(rows(cORDERED.C),1);
    
    cORDERED.CovParType 	= 1;
    cORDERED.useThreads		= 1;
    cORDERED.Algorithm 		= 4;
    cORDERED.hessMethod		= 1;
    
    if tol ne -1;
        cORDERED.DirTol = tol;
    endif;
    
    if print_it eq 1;
        cORDERED.printIters = 1;
        cORDERED.GradCheck	= 1;
    else;
        cORDERED.printIters = 0;
    endif;
    
    struct cmlmtResults ORDERED_out;
    
    if print_it eq 1;
        ORDERED_out = CMLmtprt(CMLmt(&ORDERED_MLE,pORDERED,dORDERED,cORDERED));
    else;
        ORDERED_out = CMLmt(&ORDERED_MLE,pORDERED,dORDERED,cORDERED);
    endif;
    
    beta_Y  	= pvUnpack(ORDERED_out.Par,"beta_Y");
    op_mu   = pvUnpack(ORDERED_out.Par,"mu");
    L_OP 	= ORDERED_out.Fct;
    V_OP    = real(ORDERED_out.CovPar);
    retOPb  = ORDERED_out.Retcode;
    et_OP   = ORDERED_out.ElapsedTime;
    
    OP_B  	= beta_Y|op_mu;
    
    beta_struct = beta_Y;
    
    
    if ismiss(OP_B) eq 1;
        
        OP_B  	= ones(rows(start_OP),1).*-999;
        OPbse   = OP_B;
        t_OP    = OP_B;
        L_OP	= -999;
        
    else;
        
        opbse   = sqrt(diag(V_OP));
        if rows(opbse) eq 1;
            opbse = opbse .* ones(rows(start_OP),1);
        endif;
        opbse   = missrv(opbse,10000);
        t_OP    = abs(OP_B ./ opbse);
        
    endif;
    
    OP_B = OP_B~opbse~t_OP;
    
    
    /* DETERMINISTIC PART & PROBABILITIES */
    
    zbetahat 	= X[.,2:kx_struct]*beta_Y;
    mu_mat		= op_mu' .* ones(rows(d_ij),capJ-1);
    mu_xbeta	= mu_mat - (zbetahat .*. ones(1,capJ-1));
    p_ij		= cdfn(mu_xbeta);
    p_ij		= p_ij[.,1]~(p_ij[.,2:cols(p_ij)] - p_ij[.,1:cols(p_ij)-1]);
    pordp		= p_ij~(1-sumr(p_ij));
    
    check = meanc(sumc(pordp'));
    if check gt 1.00001 or check lt (1-0.00001);
        errorlog "ORDERED PROBIT PROBABILITIES DON'T SUM TO 1: AVERAGE PROB IS";
        check;
        stop;
    endif;
    
    y_star_OP	 = zbetahat;
    
    
    /* CALL ME PROCEDURE AND CALCULATE INFORMATION CRITERIA */
    
    {IC_OP} 	= _IC(L_OP,N,rows(OP_B));
    
    
    if retopb eq 0;
        
        let retopb = normal convergence;
        
    else;
        
        let retopb 	= abnormal convergence;
        L_OP 		= -999;
        IC_OP   	= -999 .* ones(1,4);
        
    endif;
    let me_model = ordered;
    
    {me_OP, me_OPse,
        pOPxbar, pOP_se} = mes(V_OP,OP_B[.,1],capJ,me_model,meanc(X_struct[.,2:kx_struct]));
    
    if mtype eq 1;
        goto finish;
    endif;
    
    /* ESTIMATE HOPIT VIGNETTES MODEL */
    
    count_hop = 0;
    next_try:
        
        if rows(start_HOP) eq 1;
        
        beta_Y  		= OLS_B[.,1];
        
        V_cons		= rndn(kv,1) ./ 100;
        sig_V_inv	= 1;
        
        if exp_mu eq 1;
            
            gama0       = (rndn(kx_mu-1,1) ./ 100);
            gama		= (rndn(kx_mu,capJ-2) ./ 100);
            gama	    = gama';
            
        elseif exp_mu eq 0;
            
            gama = rndn(kx_mu,1) ./ 10;
            gama[1,1] = 0;
            
            for j (2,capJ-1,1);
                
                gama_temp	= gama[.,j-1] + (rndu(kx_mu,1) ./ 10);
                gama 		= gama~gama_temp;
                
            endfor;
            
            gama	    = gama';
            gama0       = gama[1,2:cols(gama)]';
            gama        = gama[2:rows(gama),.];
            
        endif;
        
    else;
        
        beta_Y  	= start_HOP[1:kx_struct];
        V_cons		= start_HOP[kx_struct+1:kx_struct+kv];
        gama0       = start_HOP[kx_struct+kv+1:kx_struct+kv+kx_mu-1];
        gama		= start_HOP[kx_struct+kv+kx_mu:kx_struct+kv+kx_mu-1+(kx_mu*(capJ-2))];
        gama		= reshape(gama,capJ-2,kx_mu);
        if sig_restrict eq 0;
            sig_V_inv	= start_HOP[kx_struct+kv+kx_mu-1+(kx_mu*(capJ-2))+1];
        endif;
        
    endif;
    
    struct PV pHOPIT;
    pHOPIT 	= pvCreate;
    pHOPIT 	= pvPack(pHOPIT,beta_Y,"beta_Y");
    pHOPIT 	= pvPack(pHOPIT,V_cons,"V_cons");
    pHOPIT 	= pvPack(pHOPIT,gama0,"gama0");
    pHOPIT 	= pvPack(pHOPIT,gama,"gama");
    if sig_restrict eq 0;
        
        pHOPIT 	= pvPack(pHOPIT,sig_V_inv,"sig_V_inv");
        k_HOPIT	= rows(vec(beta_Y)) + rows(vec(V_cons)) + rows(vec(gama0)) + rows(vec(gama)) + rows(vec(sig_V_inv));
        
    else;
        
        k_HOPIT	= rows(vec(beta_Y)) + rows(vec(V_cons)) + rows(vec(gama0)) + rows(vec(gama));
        
    endif;
    
    struct DS dHOPIT;
    if kv eq 1;
        
        dHOPIT = reshape(dsCreate,4,1);
        dHOPIT[1].DataMatrix = d_ij;
        dHOPIT[2].DataMatrix = X_struct;
        dHOPIT[3].DataMatrix = X_mu;
        dHOPIT[4].DataMatrix = V_ij1;
        
    elseif kv eq 2;
        
        dHOPIT = reshape(dsCreate,5,1);
        dHOPIT[1].DataMatrix = d_ij;
        dHOPIT[2].DataMatrix = X_struct;
        dHOPIT[3].DataMatrix = X_mu;
        dHOPIT[4].DataMatrix = V_ij1;
        dHOPIT[5].DataMatrix = V_ij2;
        
        
    elseif kv eq 3;
        
        dHOPIT = reshape(dsCreate,6,1);
        dHOPIT[1].DataMatrix = d_ij;
        dHOPIT[2].DataMatrix = X_struct;
        dHOPIT[3].DataMatrix = X_mu;
        dHOPIT[4].DataMatrix = V_ij1;
        dHOPIT[5].DataMatrix = V_ij2;
        dHOPIT[6].DataMatrix = V_ij3;
        
    endif;
    
    struct cmlmtControl cHOPIT;
    cHOPIT          = cmlmtControlCreate;
    cHOPIT.title    = "HOPIT Ordered Probit Model";
    if sig_restrict eq 0;
        cHOPIT.Bounds   = ((ones(k_HOPIT-1,1) .* -10000)~(ones(k_HOPIT-1,1) .* 10000))|(0.00001~10);
    endif;
    cHOPIT.useThreads	= 1;
    if kv eq 1;
        CHOPIT.MaxIters		= 1000;
    elseif kv eq 2;
        CHOPIT.MaxIters		= 2000;
    elseif kv eq 3;
        CHOPIT.MaxIters		= 3000;
    endif;
    
    if tol ne -1;
        cHOPIT.DirTol = tol;
    endif;
    
    if print_it eq 1;
        
        cHOPIT.printIters 	= 1;
        cHOPIT.GradCheck	= 1;
        
    else;
        
        cHOPIT.printIters = 0;
        
    endif;
    
    struct cmlmtResults HOPIT_out;
    
    HOPIT_fail = 0;
    
    if print_it eq 1;
        HOPIT_out = CMLmtprt(CMLmt(&HOPIT_MLE,pHOPIT,dHOPIT,cHOPIT));
    else;
        HOPIT_out = CMLmt(&HOPIT_MLE,pHOPIT,dHOPIT,cHOPIT);
    endif;
    
    //print pvGetParnames(HOPIT_out.par);
    
    beta_Y  	= pvUnpack(HOPIT_out.Par,"beta_Y");
    V_cons  	= pvUnpack(HOPIT_out.Par,"V_cons");
    gama  		= pvUnpack(HOPIT_out.Par,"gama");
    gama0  		= pvUnpack(HOPIT_out.Par,"gama0");
    if cols(gama0) eq 1;
        gama0 = gama0';
    endif;
    
    beta_HOP	= beta_Y;
    gama0_HOP	= gama0;
    gama_HOP	= gama;
    gama        = (0~gama0)|gama;
    gama		= gama';
    
    if sig_restrict eq 0;
        sig_V_inv  	= pvUnpack(HOPIT_out.Par,"sig_V_inv");
    else;
        sig_V_inv = 1;
    endif;
    
    L_HOPIT 	= HOPIT_out.Fct;
    V_HOPIT    	= real(HOPIT_out.CovPar);
    retHOPIT  	= HOPIT_out.Retcode;
    et_HOPIT   	= HOPIT_out.ElapsedTime;
    HOPIT_grad  = HOPIT_out.Gradient;
    
    if ismiss(V_HOPIT);
        
        HOPIT_fail  = 1;
        
        if det(HOPIT_out.xproduct) gt 0;
            
            V_HOPIT = eye(rows(HOPIT_out.xproduct))/HOPIT_out.xproduct;
            
        else;
            
            V_HOPIT     = 1;
            V_HOPIT     = miss(V_HOPIT,1);
            HOPIT_fail  = 1;
            
        endif;
        
    endif;
    
    if sig_V_inv ne 1;
        
        if ismiss(V_HOPIT) ne 1;
            V_gamma     = V_HOPIT[rows(beta_Y)+rows(V_cons)+1:rows(beta_Y)+rows(V_cons)+rows(vec(gama)),rows(beta_Y)+rows(V_cons)+1:rows(beta_Y)+rows(V_cons)+rows(vec(gama))];
        endif;
        
    else;
        
        if ismiss(V_HOPIT) ne 1;
            V_gamma     = V_HOPIT[rows(beta_Y)+rows(V_cons)+1:rows(beta_Y)+rows(V_cons)+rows(vec(gama))-1,rows(beta_Y)+rows(V_cons)+1:rows(beta_Y)+rows(V_cons)+rows(vec(gama))-1];
        endif;
        
    endif;
    
    if sig_restrict eq 0;
        HOPIT_B   	= beta_Y|V_cons|vec(gama)|sig_V_inv;
    else;
        HOPIT_B   	= beta_Y|V_cons|vec(gama);
    endif;
    
    if maxc(abs(HOPIT_B)) ge 10;
        HOPIT_fail = 1;
    endif;
    
    if ismiss(HOPIT_B) eq 1;
        
        fail:
            if fail_linear eq 1;
            stop;
        endif;
        HOPIT_B		= ones(rows(start_HOP),1) .* 10000;
        HOPIT_se	= HOPIT_B;
        t_HOPIT    	= HOPIT_B;
        L_HOPIT		= -999;
        LR_test     = -999~-999~-999;
        LM_test     = (-999~-999~-999)|(-999~-999~-999)|(-999~-999~-999);
        Wald_test	= -999~-999~-999;
        HOPIT_fail  = 1;
        
    else;
        
        HOPIT_se   	= sqrt(diag(V_HOPIT));
        
        if rows(HOPIT_se) eq 1;
            
            HOPIT_se    = HOPIT_se .* ones(rows(HOPIT_B),1);
            HOPIT_fail  = 1;
            
        endif;
        
        HOPIT_se   	= missrv(HOPIT_se,1000);
        
        if sumc(HOPIT_se .eq 10000) ge 1;
            
            HOPIT_fail  = 1;
            
        endif;
        
        if rows(HOPIT_se) ne k_HOPIT;
            
            HOPIT_fail  = 1;
            
        endif;
        
        if rows(V_HOPIT) eq 1;
            
            HOPIT_se    = HOPIT_se .* ones(k_HOPIT,1);
            HOPIT_fail  = 1;
            
        endif;
        
        if sig_restrict eq 0;
            se_gama		= -999|HOPIT_se[rows(beta_Y)+rows(V_cons)+1:rows(HOPIT_se)-1];
        elseif sig_restrict eq 1;
            se_gama		= -999|HOPIT_se[rows(beta_Y)+rows(V_cons)+1:rows(HOPIT_se)];
        endif;
        
        se_gama		= reshape(se_gama,kx_mu,capJ-1)';
        se_gama		= vec(se_gama);
        
        if sig_V_inv ne 1;
            HOPIT_se 	= HOPIT_se[1:rows(beta_Y)+rows(V_cons)]|se_gama|HOPIT_se[rows(beta_Y)+rows(V_cons)+rows(se_gama):rows(HOPIT_se)];
        else;
            HOPIT_se 	= HOPIT_se[1:rows(beta_Y)+rows(V_cons)]|se_gama;
        endif;
        
        t_HOPIT    	= abs(HOPIT_B ./ HOPIT_se);
        
    endif;
    
    HOPIT_B   		= HOPIT_B~HOPIT_se~abs(HOPIT_B ./ HOPIT_se);
    
    /* DEFINE COMMON BOUNDARY PARAMETERS */
    
    gama_LR    = gama;
    beta_LR    = beta_Y;
    zgama	    = X_mu*gama;
    if exp_mu eq 1;
        mu_ij	= _mu_ij(zgama);
    elseif exp_mu eq 0;
        mu_ij	= zgama;
    endif;
    
    /* DETERMINISTIC PART & PROBABILITIES */
    
    xbeta 			= X_struct*beta_Y;
    mu_xbeta		= mu_ij - (xbeta .*. ones(1,capJ-1));
    p_ij			= cdfn(mu_xbeta);
    p_ij			= p_ij[.,1]~(p_ij[.,2:cols(p_ij)] - p_ij[.,1:cols(p_ij)-1]);
    p_HOPIT			= p_ij~(1-sumr(p_ij));
    
    mu0_star = mu_xbeta[.,1];
    mu1_star = mu0_star + mu_ij[.,2];
    if capJ gt 3;
        mu2_star = mu1_star + mu_ij[.,3];
    endif;
    
    /* check for ill-defined probs */
    
    if minc(vec(p_HOPIT)) lt 0;
        
        HOPIT_fail  = 1;
        fail_linear = 1;
        
    endif;
    
    mu_xbetaV1	= (mu_ij - V_cons[1]) .* sig_V_inv;
    p_ijV1		= cdfn(mu_xbetaV1);
    p_ijV1		= p_ijV1[.,1]~(p_ijV1[.,2:cols(p_ijV1)] - p_ijV1[.,1:cols(p_ijV1)-1]);
    p_ijV1		= p_ijV1~(1-sumr(p_ijV1));
    
    if minc(vec(p_ijV1)) lt 0;
        
        HOPIT_fail  = 1;
        fail_linear = 1;
        
    endif;
    
    if kv ge 2;
        
        mu_xbetaV2	= (mu_ij - V_cons[2]) .* sig_V_inv;
        p_ijV2		= cdfn(mu_xbetaV2);
        p_ijV2		= p_ijV2[.,1]~(p_ijV2[.,2:cols(p_ijV2)] - p_ijV2[.,1:cols(p_ijV2)-1]);
        p_ijV2		= p_ijV2~(1-sumr(p_ijV2));
        
        if minc(vec(p_ijV2)) lt 0;
            
            HOPIT_fail  = 1;
            fail_linear = 1;
            
        endif;
        
    endif;
    
    if kv ge 3;
        
        mu_xbetaV3	= (mu_ij - V_cons[3]) .* sig_V_inv;
        p_ijV3		= cdfn(mu_xbetaV3);
        p_ijV3		= p_ijV3[.,1]~(p_ijV3[.,2:cols(p_ijV3)] - p_ijV3[.,1:cols(p_ijV3)-1]);
        p_ijV3		= p_ijV3~(1-sumr(p_ijV3));
        
        if minc(vec(p_ijV3)) lt 0;
            
            HOPIT_fail  = 1;
            fail_linear = 1;
            
        endif;
        
    endif;
    
    lnLiHOPIT		= ln(sumr(d_ij .* p_HOPIT));
    
    y_star_HOPIT	= xbeta;
    
    /* CALL ME PROCEDURE AND CALCULATE INFORMATION CRITERIA */
    
    {IC_HOPIT} 	= _IC(L_HOPIT,N,k_HOPIT);
    
    if retHOPIT eq 0;
        
        let retHOPIT = normal convergence;
        
    else;
        
        let retHOPIT 	= abnormal convergence;
        HOPIT_fail      = 1;
        
    endif;
    
    //HOPIT_fail;stop;
    
    let me_model 	= hopit;
    
    {me_HOPIT, me_HOPITse,
        p_HOPIT_xbar, p_HOPITse} = mes(V_HOPIT,delif(HOPIT_B[.,1], HOPIT_B[.,1] .eq 0),capJ,me_model,meanc(X));
    
    /* VIGNETTE SPECIFICATION TESTS */
    
    if exp_mu eq 0;
        
    elseif exp_mu eq 1;
        
        if LR_Yes_No ge 1;
            
            if HOPIT_fail ge 1;
                
                LM_test     = (-999~-999~-999)|(-999~-999~-999)|(-999~-999~-999);
                goto skip_LM;
                
            endif;
            
            /* LM TEST */
            
            ?;
            "max abs average gradient";;
            maxc(meanc(abs(HOPIT_grad)));
            ?;
            
            if maxc(meanc(abs(HOPIT_grad))) ge 0.01;
                
                LM_test = ones(3,3) .* -999;
                goto skip_LM;
                
            endif;
            
            /* define derivatives under the null */
            
            beta_V1  	= V_cons[1]|zeros(kx_struct-1,1);
            if kv ge 2;
                beta_V2     = V_cons[2]|zeros(kx_struct-1,1);
                if kv ge 3;
                    beta_V3     = V_cons[3]|zeros(kx_struct-1,1);
                endif;
            endif;
            
            gama_V1		= gama;
            if kv ge 2;
                gama_V2 = gama;
                if kv ge 3;
                    gama_V3 = gama;
                endif;
            endif;
            
            xbeta_V1 	= X_struct*beta_V1;
            mu_xbeta_V1	= (mu_ij - (xbeta_V1 .*. ones(1,capJ-1))) .* sig_V_inv;
            p_ij_V1		= cdfn(mu_xbeta_V1);
            p_ij_V1		= p_ij_V1[.,1]~(p_ij_V1[.,2:cols(p_ij_V1)] - p_ij_V1[.,1:cols(p_ij_V1)-1]);
            p_ij_V1		= p_ij_V1~(1-sumr(p_ij_V1));
            
            if kv ge 2;
                
                xbeta_V2 	= X_struct*beta_V2;
                mu_xbeta_V2	= (mu_ij - (xbeta_V2 .*. ones(1,capJ-1))) .* sig_V_inv;
                p_ij_V2		= cdfn(mu_xbeta_V2);
                p_ij_V2		= p_ij_V2[.,1]~(p_ij_V2[.,2:cols(p_ij_V2)] - p_ij_V2[.,1:cols(p_ij_V2)-1]);
                p_ij_V2		= p_ij_V2~(1-sumr(p_ij_V2));
                
                if kv ge 3;
                    
                    xbeta_V3 	= X_struct*beta_V3;
                    mu_xbeta_V3	= (mu_ij - (xbeta_V3 .*. ones(1,capJ-1))) .* sig_V_inv;
                    p_ij_V3		= cdfn(mu_xbeta_V3);
                    p_ij_V3		= p_ij_V3[.,1]~(p_ij_V3[.,2:cols(p_ij_V3)] - p_ij_V3[.,1:cols(p_ij_V3)-1]);
                    p_ij_V3		= p_ij_V3~(1-sumr(p_ij_V3));
                    
                endif;
                
            endif;
            
            /* theta = 1/sigma derivatives; V1, V2 and V3 */
            
            if sig_restrict eq 0;
                
                if kv eq 1;
                    
                    {grad_theta_V1}	= _grad_theta_HOPIT(mu_xbetaV1,sig_V_inv,p_ijV1,V_ij1);
                    grad_theta		= grad_theta_V1;
                    
                elseif kv eq 2;
                    
                    {grad_theta_V1}	= _grad_theta_HOPIT(mu_xbetaV1,sig_V_inv,p_ijV1,V_ij1);
                    {grad_theta_V2}	= _grad_theta_HOPIT(mu_xbetaV2,sig_V_inv,p_ijV2,V_ij2);
                    grad_theta		= grad_theta_V1 + grad_theta_V2;
                    
                elseif kv eq 3;
                    
                    {grad_theta_V1}	= _grad_theta_HOPIT(mu_xbetaV1,sig_V_inv,p_ijV1,V_ij1);
                    {grad_theta_V2}	= _grad_theta_HOPIT(mu_xbetaV2,sig_V_inv,p_ijV2,V_ij2);
                    {grad_theta_V3}	= _grad_theta_HOPIT(mu_xbetaV3,sig_V_inv,p_ijV3,V_ij3);
                    grad_theta		= grad_theta_V1 + grad_theta_V2 + grad_theta_V3;
                    
                endif;
                
            endif;
            
            
            if capJ eq 5;
                
                /* beta_Y derivatives */
                {grad_beta} = _grad_beta_OP(mu_xbeta,d_ij,p_HOPIT,X_struct);
                
                /* gradient for gamma_1 HOPIT */
                {grad_g1}	= _grad_g_j(gama[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta,p_HOPIT,d_ij,1);
                
                /* gradient for gamma_2 HOPIT */
                {grad_g2}	= _grad_g_j(gama[.,2],X_mu,mu_xbeta[.,2:cols(mu_xbeta)],p_HOPIT[.,2:capJ],d_ij[.,2:capJ],1);
                
                /* gradient for gamma_3 HOPIT */
                {grad_g3}	= _grad_g_j(gama[.,3],X_mu,mu_xbeta[.,3:cols(mu_xbeta)],p_HOPIT[.,3:capJ],d_ij[.,3:capJ],1);
                
                /* gradient for gamma_4 HOPIT */
                {grad_g4}	= _grad_gJ(gama[.,4],X_mu,mu_xbeta[.,cols(mu_xbeta)],p_HOPIT[.,4:capJ],d_ij[.,4:capJ],1);
                
                //Grad 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_g4;
                
            elseif capJ eq 4;
                
                /* beta_Y derivatives */
                {grad_beta} = _grad_beta_OP(mu_xbeta,d_ij,p_HOPIT,X_struct);
                
                /* gradient for gamma_1 HOPIT */
                {grad_g1}	= _grad_g_j(gama[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta,p_HOPIT,d_ij,1);
                
                /* gradient for gamma_2 HOPIT */
                {grad_g2}	= _grad_g_j(gama[.,2],X_mu,mu_xbeta[.,2:cols(mu_xbeta)],p_HOPIT[.,2:capJ],d_ij[.,2:capJ],1);
                
                /* gradient for gamma_3 HOPIT */
                {grad_g3}	= _grad_gJ(gama[.,3],X_mu,mu_xbeta[.,cols(mu_xbeta)],p_HOPIT[.,3:capJ],d_ij[.,3:capJ],1);
                
                //Grad 	= grad_beta~grad_g1~grad_g2~grad_g3;
                
            elseif capJ eq 3;
                
                /* beta_Y derivatives */
                {grad_beta} = _grad_beta_OP(mu_xbeta,d_ij,p_HOPIT,X_struct);
                
                /* gradient for gamma_1 HOPIT */
                {grad_g1}	= _grad_g_j(gama[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta,p_HOPIT,d_ij,1);
                
                /* gradient for gamma_3 HOPIT */
                {grad_g2}	= _grad_gJ(gama[.,2],X_mu,mu_xbeta[.,cols(mu_xbeta)],p_HOPIT[.,2:capJ],d_ij[.,2:capJ],1);
                
            endif;
            
            if capJ eq 5;
                
                /* beta_V1 derivatives */
                {grad_beta_V1} = _grad_beta_OP(mu_xbeta_V1,V_ij1,p_ij_V1,X_struct);
                grad_beta_V1 	= grad_beta_V1 .* sig_V_inv;
                
                /* gradient for gamma_1_V1 HOPIT */
                {grad_g1_V1}	= _grad_g_j(gama_V1[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta_V1,p_ij_V1,V_ij1,1);
                grad_g1_V1		= grad_g1_V1 .* sig_V_inv;
                
                /* gradient for gamma_2_V1 HOPIT */
                {grad_g2_V1}	= _grad_g_j(gama_V1[.,2],X_mu,mu_xbeta_V1[.,2:cols(mu_xbeta_V1)],p_ij_V1[.,2:capJ],V_ij1[.,2:capJ],1);
                grad_g2_V1		= grad_g2_V1 .* sig_V_inv;
                
                /* gradient for gamma_3_V1 HOPIT */
                
                {grad_g3_V1}	= _grad_g_j(gama_V1[.,3],X_mu,mu_xbeta_V1[.,3:cols(mu_xbeta_V1)],p_ij_V1[.,3:capJ],V_ij1[.,3:capJ],1);
                grad_g3_V1		= grad_g3_V1 .* sig_V_inv;
                
                /* gradient for gamma_4_V1 HOPIT */
                {grad_g4_V1}	= _grad_gJ(gama_V1[.,4],X_mu,mu_xbeta_V1[.,cols(mu_xbeta_V1)],p_ij_V1[.,4:capJ],V_ij1[.,4:capJ],1);
                grad_g4_V1		= grad_g4_V1 .* sig_V_inv;
                
                if kv ge 2;
                    
                    /* beta_V2 derivatives */
                    {grad_beta_V2} = _grad_beta_OP(mu_xbeta_V2,V_ij2,p_ij_V2,X_struct);
                    grad_beta_V2 	= grad_beta_V2 .* sig_V_inv;
                    
                    /* gradient for gamma_1_V2 HOPIT */
                    {grad_g1_V2}	= _grad_g_j(gama_V2[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta_V2,p_ij_V2,V_ij2,1);
                    grad_g1_V2		= grad_g1_V2 .* sig_V_inv;
                    
                    /* gradient for gamma_2_V2 HOPIT */
                    {grad_g2_V2}	= _grad_g_j(gama_V2[.,2],X_mu,mu_xbeta_V2[.,2:cols(mu_xbeta_V2)],p_ij_V2[.,2:capJ],V_ij2[.,2:capJ],1);
                    grad_g2_V2		= grad_g2_V2 .* sig_V_inv;
                    
                    /* gradient for gamma_3_V2 HOPIT */
                    
                    {grad_g3_V2}	= _grad_g_j(gama_V2[.,3],X_mu,mu_xbeta_V2[.,3:cols(mu_xbeta_V2)],p_ij_V2[.,3:capJ],V_ij2[.,3:capJ],1);
                    grad_g3_V2		= grad_g3_V2 .* sig_V_inv;
                    
                    /* gradient for gamma_4_V2 HOPIT */
                    {grad_g4_V2}	= _grad_gJ(gama_V2[.,4],X_mu,mu_xbeta_V2[.,cols(mu_xbeta_V2)],p_ij_V2[.,4:capJ],V_ij2[.,4:capJ],1);
                    grad_g4_V2		= grad_g4_V2 .* sig_V_inv;
                    
                    if kv ge 3;
                        
                        /* beta_V3 derivatives */
                        {grad_beta_V3} = _grad_beta_OP(mu_xbeta_V3,V_ij3,p_ij_V3,X_struct);
                        grad_beta_V3 	= grad_beta_V3 .* sig_V_inv;
                        
                        /* gradient for gamma_1_V3 HOPIT */
                        {grad_g1_V3}	= _grad_g_j(gama_V3[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta_V3,p_ij_V3,V_ij3,1);
                        grad_g1_V3		= grad_g1_V3 .* sig_V_inv;
                        
                        /* gradient for gamma_2_V3 HOPIT */
                        {grad_g2_V3}	= _grad_g_j(gama_V3[.,2],X_mu,mu_xbeta_V3[.,2:cols(mu_xbeta_V3)],p_ij_V3[.,2:capJ],V_ij3[.,2:capJ],1);
                        grad_g2_V3		= grad_g2_V3 .* sig_V_inv;
                        
                        /* gradient for gamma_3_V3 HOPIT */
                        
                        {grad_g3_V3}	= _grad_g_j(gama_V3[.,3],X_mu,mu_xbeta_V3[.,3:cols(mu_xbeta_V3)],p_ij_V3[.,3:capJ],V_ij3[.,3:capJ],1);
                        grad_g3_V3		= grad_g3_V3 .* sig_V_inv;
                        
                        /* gradient for gamma_4_V3 HOPIT */
                        {grad_g4_V3}	= _grad_gJ(gama_V3[.,4],X_mu,mu_xbeta_V3[.,cols(mu_xbeta_V3)],p_ij_V3[.,4:capJ],V_ij3[.,4:capJ],1);
                        grad_g4_V3		= grad_g4_V3 .* sig_V_inv;
                        
                    endif;
                    
                endif;
                
                if kv eq 1;
                    Score_anal 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_g4~grad_beta_V1~grad_g1_V1~grad_g2_V1~grad_g3_V1~grad_g4_V1;
                elseif kv eq 2;
                    Score_anal 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_g4~grad_beta_V1~grad_g1_V1~grad_g2_V1~grad_g3_V1~grad_g4_V1~grad_beta_V2~grad_g1_V2~grad_g2_V2~grad_g3_V2~grad_g4_V2;
                elseif kv eq 3;
                    Score_anal 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_g4~grad_beta_V1~grad_g1_V1~grad_g2_V1~grad_g3_V1~grad_g4_V1~grad_beta_V2~grad_g1_V2~grad_g2_V2~grad_g3_V2~grad_g4_V2~
                        grad_beta_V3~grad_g1_V3~grad_g2_V3~grad_g3_V3~grad_g4_V3;
                endif;
                
                if kv eq 1;
                    Score_anal_VE 	= grad_beta~(grad_g1+grad_g1_V1)~(grad_g2+grad_g2_V1)~(grad_g3+grad_g3_V1)~(grad_g4+grad_g4_V1)~grad_beta_V1;
                elseif kv eq 2;
                    Score_anal_VE 	= grad_beta~(grad_g1+grad_g1_V1+grad_g1_V2)~(grad_g2+grad_g2_V1+grad_g2_V2)~(grad_g3+grad_g3_V1+grad_g3_V2)~(grad_g4+grad_g4_V1+grad_g4_V2)~grad_beta_V1~grad_beta_V2;
                elseif kv eq 3;
                    Score_anal_VE 	= grad_beta~(grad_g1+grad_g1_V1+grad_g1_V2+grad_g1_V3)~(grad_g2+grad_g2_V1+grad_g2_V2+grad_g2_V3)~(grad_g3+grad_g3_V1+grad_g3_V2+grad_g3_V3)~
                        (grad_g4+grad_g4_V1+grad_g4_V2+grad_g4_V3)~grad_beta_V1~grad_beta_V2~grad_beta_V2;
                endif;
                
                if kv eq 1;
                    Score_anal_RC 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_g4~grad_beta_V1[.,1]~grad_g1_V1~grad_g2_V1~grad_g3_V1~grad_g4_V1;
                elseif kv eq 2;
                    Score_anal_RC 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_g4~grad_beta_V1[.,1]~grad_g1_V1[.,1]~grad_g2_V1~grad_g3_V1~grad_g4_V1~grad_beta_V2~grad_g1_V2~grad_g2_V2~grad_g3_V2~grad_g4_V2;
                elseif kv eq 3;
                    Score_anal_RC 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_g4~grad_beta_V1[.,1]~grad_g1_V1~grad_g2_V1~grad_g3_V1~grad_g4_V1~grad_beta_V2[.,1]~grad_g1_V2~grad_g2_V2~grad_g3_V2~grad_g4_V2~
                        grad_beta_V3[.,1]~grad_g1_V3~grad_g2_V3~grad_g3_V3~grad_g4_V3;
                endif;
                
            elseif capJ eq 4;
                
                /* beta_V1 derivatives */
                {grad_beta_V1} = _grad_beta_OP(mu_xbeta_V1, V_ij1, p_ij_V1, X_struct);
                grad_beta_V1 	= grad_beta_V1 .* sig_V_inv;
                
                /* gradient for gamma_1 HOPIT */
                {grad_g1_V1}	= _grad_g_j(gama_V1[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta_V1,p_ij_V1,V_ij1,1);
                grad_g1_V1		= grad_g1_V1 .* sig_V_inv;
                
                /* gradient for gamma_2 HOPIT */
                {grad_g2_V1}	= _grad_g_j(gama_V1[.,2],X_mu,mu_xbeta_V1[.,2:cols(mu_xbeta_V1)],p_ij_V1[.,2:capJ],V_ij1[.,2:capJ],1);
                grad_g2_V1		= grad_g2_V1 .* sig_V_inv;
                
                /* gradient for gamma_3 HOPIT */
                {grad_g3_V1}	= _grad_gJ(gama_V1[.,3],X_mu,mu_xbeta_V1[.,cols(mu_xbeta_V1)],p_ij_V1[.,3:capJ],V_ij1[.,3:capJ],1);
                grad_g3_V1		= grad_g3_V1 .* sig_V_inv;
                
                if kv ge 2;
                    
                    /* beta_V2 derivatives */
                    {grad_beta_V2} = _grad_beta_OP(mu_xbeta_V2,V_ij2,p_ij_V2,X_struct);
                    grad_beta_V2 	= grad_beta_V2 .* sig_V_inv;
                    
                    /* gradient for gamma_1 HOPIT */
                    {grad_g1_V2}	= _grad_g_j(gama_V2[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta_V2,p_ij_V2,V_ij2,1);
                    grad_g1_V2		= grad_g1_V2 .* sig_V_inv;
                    
                    /* gradient for gamma_2 HOPIT */
                    {grad_g2_V2}	= _grad_g_j(gama_V2[.,2],X_mu,mu_xbeta_V2[.,2:cols(mu_xbeta_V2)],p_ij_V2[.,2:capJ],V_ij2[.,2:capJ],1);
                    grad_g2_V2		= grad_g2_V2 .* sig_V_inv;
                    
                    /* gradient for gamma_3 HOPIT */
                    {grad_g3_V2}	= _grad_gJ(gama_V2[.,3],X_mu,mu_xbeta_V2[.,cols(mu_xbeta_V2)],p_ij_V2[.,3:capJ],V_ij2[.,3:capJ],1);
                    grad_g3_V2		= grad_g3_V2 .* sig_V_inv;
                    
                    if kv ge 3;
                        
                        /* beta_V3 derivatives */
                        {grad_beta_V3} = _grad_beta_OP(mu_xbeta_V3,V_ij3,p_ij_V3,X_struct);
                        grad_beta_V3 	= grad_beta_V3 .* sig_V_inv;
                        
                        /* gradient for gamma_1 HOPIT */
                        {grad_g1_V3}	= _grad_g_j(gama_V3[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta_V3,p_ij_V3,V_ij3,1);
                        grad_g1_V3		= grad_g1_V3 .* sig_V_inv;
                        
                        /* gradient for gamma_2 HOPIT */
                        {grad_g2_V3}	= _grad_g_j(gama_V3[.,2],X_mu,mu_xbeta_V3[.,2:cols(mu_xbeta_V3)],p_ij_V3[.,2:capJ],V_ij3[.,2:capJ],1);
                        grad_g2_V3		= grad_g2_V3 .* sig_V_inv;
                        
                        /* gradient for gamma_3 HOPIT */
                        {grad_g3_V3}	= _grad_gJ(gama_V3[.,3],X_mu,mu_xbeta_V3[.,cols(mu_xbeta_V3)],p_ij_V3[.,3:capJ],V_ij3[.,3:capJ],1);
                        grad_g3_V3		= grad_g3_V3 .* sig_V_inv;
                        
                    endif;
                    
                endif;
                
                if kv eq 1;
                    Score_anal 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_beta_V1~grad_g1_V1~grad_g2_V1~grad_g3_V1;
                elseif kv eq 2;
                    Score_anal 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_beta_V1~grad_g1_V1~grad_g2_V1~grad_g3_V1~grad_beta_V2~grad_g1_V2~grad_g2_V2~grad_g3_V2;
                elseif kv eq 3;
                    Score_anal 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_beta_V1~grad_g1_V1~grad_g2_V1~grad_g3_V1~grad_beta_V2~grad_g1_V2~grad_g2_V2~grad_g3_V2~
                        grad_beta_V3~grad_g1_V3~grad_g2_V3~grad_g3_V3;
                endif;
                
                if kv eq 1;
                    Score_anal_VE 	= grad_beta~(grad_g1+grad_g1_V1)~(grad_g2+grad_g2_V1)~(grad_g3+grad_g3_V1)~grad_beta_V1;
                elseif kv eq 2;
                    Score_anal_VE 	= grad_beta~(grad_g1+grad_g1_V1+grad_g1_V2)~(grad_g2+grad_g2_V1+grad_g2_V2)~(grad_g3+grad_g3_V1+grad_g3_V2)~grad_beta_V1~grad_beta_V2;
                elseif kv eq 3;
                    Score_anal_VE 	= grad_beta~(grad_g1+grad_g1_V1+grad_g1_V2+grad_g1_V3)~(grad_g2+grad_g2_V1+grad_g2_V2+grad_g2_V3)~(grad_g3+grad_g3_V1+grad_g3_V2+grad_g3_V3)~
                        grad_beta_V1~grad_beta_V2~grad_beta_V3;
                endif;
                
                if kv eq 1;
                    Score_anal_RC 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_beta_V1[.,1]~grad_g1_V1~grad_g2_V1~grad_g3_V1;
                elseif kv eq 2;
                    Score_anal_RC 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_beta_V1[.,1]~grad_g1_V1~grad_g2_V1~grad_g3_V1~grad_beta_V2[.,1]~grad_g1_V2~grad_g2_V2~grad_g3_V2;
                elseif kv eq 3;
                    Score_anal_RC 	= grad_beta~grad_g1~grad_g2~grad_g3~grad_beta_V1[.,1]~grad_g1_V1~grad_g2_V1~grad_g3_V1~grad_beta_V2[.,1]~grad_g1_V2~grad_g2_V2~grad_g3_V2~
                        grad_beta_V3[.,1]~grad_g1_V3~grad_g2_V3~grad_g3_V3;
                endif;
                
                //END OF J=4 BLOCK
                
            elseif capJ eq 3;
                
                /* beta_V1 derivatives */
                {grad_beta_V1} = _grad_beta_OP(mu_xbeta_V1,V_ij1,p_ij_V1,X_struct);
                grad_beta_V1 	= grad_beta_V1 .* sig_V_inv;
                
                /* gradient for gamma_1 HOPIT */
                {grad_g1_V1}	= _grad_g_j(gama_V1[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta_V1,p_ij_V1,V_ij1,1);
                grad_g1_V1		= grad_g1_V1 .* sig_V_inv;
                
                /* gradient for gamma_2 HOPIT */
                {grad_g2_V1}	= _grad_gJ(gama_V1[.,2],X_mu,mu_xbeta_V1[.,cols(mu_xbeta_V1)],p_ij_V1[.,2:capJ],V_ij1[.,2:capJ],1);
                grad_g2_V1		= grad_g2_V1 .* sig_V_inv;
                
                if kv ge 2;
                    
                    /* beta_V2 derivatives */
                    {grad_beta_V2} = _grad_beta_OP(mu_xbeta_V2,V_ij2,p_ij_V2,X_struct);
                    grad_beta_V2 	= grad_beta_V2 .* sig_V_inv;
                    
                    /* gradient for gamma_1 HOPIT */
                    {grad_g1_V2}	= _grad_g_j(gama_V2[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta_V2,p_ij_V2,V_ij2,1);
                    grad_g1_V2		= grad_g1_V2 .* sig_V_inv;
                    
                    /* gradient for gamma_2 HOPIT */
                    {grad_g2_V2}	= _grad_gJ(gama_V2[.,2],X_mu,mu_xbeta_V2[.,cols(mu_xbeta_V2)],p_ij_V2[.,2:capJ],V_ij2[.,2:capJ],1);
                    grad_g2_V2		= grad_g2_V2 .* sig_V_inv;
                    
                    if kv ge 3;
                        
                        /* beta_V3 derivatives */
                        {grad_beta_V3} = _grad_beta_OP(mu_xbeta_V3,V_ij3,p_ij_V3,X_struct);
                        grad_beta_V3 	= grad_beta_V3 .* sig_V_inv;
                        
                        /* gradient for gamma_1 HOPIT */
                        {grad_g1_V3}	= _grad_g_j(gama_V3[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta_V3,p_ij_V3,V_ij3,1);
                        grad_g1_V3		= grad_g1_V3 .* sig_V_inv;
                        
                        /* gradient for gamma_3 HOPIT */
                        {grad_g2_V3}	= _grad_gJ(gama_V3[.,2],X_mu,mu_xbeta_V3[.,cols(mu_xbeta_V3)],p_ij_V3[.,2:capJ],V_ij3[.,2:capJ],1);
                        grad_g2_V3		= grad_g2_V3 .* sig_V_inv;
                        
                    endif;
                    
                endif;
                
                if kv eq 1;
                    Score_anal 	= grad_beta~grad_g1~grad_g2~grad_beta_V1~grad_g1_V1~grad_g2_V1;
                elseif kv eq 2;
                    Score_anal 	= grad_beta~grad_g1~grad_g2~grad_beta_V1~grad_g1_V1~grad_g2_V1~grad_beta_V2~grad_g1_V2~grad_g2_V2;
                elseif kv eq 3;
                    Score_anal 	= grad_beta~grad_g1~grad_g2~grad_beta_V1~grad_g1_V1~grad_g2_V1~grad_beta_V2~grad_g1_V2~grad_g2_V2~
                        grad_beta_V3~grad_g1_V3~grad_g2_V3;
                endif;
                
                if kv eq 1;
                    Score_anal_VE 	= grad_beta~(grad_g1+grad_g1_V1)~(grad_g2+grad_g2_V1)~grad_beta_V1;
                elseif kv eq 2;
                    Score_anal_VE 	= grad_beta~(grad_g1+grad_g1_V1+grad_g1_V2)~(grad_g2+grad_g2_V1+grad_g2_V2)~grad_beta_V1~grad_beta_V2;
                elseif kv eq 3;
                    Score_anal_VE 	= grad_beta~(grad_g1+grad_g1_V1+grad_g1_V2+grad_g1_V3)~(grad_g2+grad_g2_V1+grad_g2_V2+grad_g2_V3)~
                        grad_beta_V1~grad_beta_V2~grad_beta_V3;
                endif;
                
                if kv eq 1;
                    Score_anal_RC 	= grad_beta~grad_g1~grad_g2~grad_beta_V1[.,1]~grad_g1_V1~grad_g2_V1;
                elseif kv eq 2;
                    Score_anal_RC 	= grad_beta~grad_g1~grad_g2~grad_beta_V1[.,1]~grad_g1_V1~grad_g2_V1~grad_beta_V2[.,1]~grad_g1_V2~grad_g2_V2;
                elseif kv eq 3;
                    Score_anal_RC 	= grad_beta~grad_g1~grad_g2~grad_beta_V1[.,1]~grad_g1_V1~grad_g2_V1~grad_beta_V2[.,1]~grad_g1_V2~grad_g2_V2~
                        grad_beta_V3[.,1]~grad_g1_V3~grad_g2_V3;
                endif;
                
                //END OF J=3 BLOCK
                
            endif;
            
            k_GOP_JOINT	= (rows(vec(beta_HOP)) + rows(vec(gama0_HOP)) + rows(vec(gama_HOP)))*(kv+1);
            
            if sig_restrict eq 0;
                
                Score_anal      = Score_anal~grad_theta;
                Score_anal_RC   = Score_anal_RC~grad_theta;
                Score_anal_VE   = Score_anal_VE~grad_theta;
                
            endif;
            
            df				= cols(Score_anal) - (k_HOPIT);
            bhhh_anal		= Score_anal'Score_anal;
            
            df_RC			= cols(score_anal_RC) - (k_HOPIT);
            bhhh_anal_RC	= score_anal_RC'score_anal_RC;
            
            df_VE			= cols(score_anal_VE) - (k_HOPIT);
            bhhh_anal_VE	= score_anal_VE'score_anal_VE;
            
            bhhh_max    	= 1;
            bhhh_anal   	= bhhh_anal ./ bhhh_max;
            bhhh_anal_RC   	= bhhh_anal_RC ./ bhhh_max;
            bhhh_anal_VE   	= bhhh_anal_VE ./ bhhh_max;
            
            if capJ eq 3;
                
                if kv eq 1;
                    score_anal 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_beta_V1)|sumc(grad_g1_V1)|sumc(grad_g2_V1));
                elseif kv eq 2;
                    score_anal 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_beta_V1)|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_beta_V2)|sumc(grad_g1_V2)|sumc(grad_g2_V2));
                elseif kv eq 3;
                    score_anal 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_beta_V1)|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_beta_V2)|sumc(grad_g1_V2)|sumc(grad_g2_V2)|
                        sumc(grad_beta_V3)|sumc(grad_g1_V3)|sumc(grad_g2_V3));
                endif;
                
            elseif capJ eq 4;
                
                if kv eq 1;
                    score_anal 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_beta_V1)|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1));
                elseif kv eq 2;
                    score_anal 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_beta_V1)|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1)|sumc(grad_beta_V2)|sumc(grad_g1_V2)|sumc(grad_g2_V2)|sumc(grad_g3_V2));
                elseif kv eq 3;
                    score_anal 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_beta_V1)|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1)|sumc(grad_beta_V2)|sumc(grad_g1_V2)|sumc(grad_g2_V2)|sumc(grad_g3_V2)|
                        sumc(grad_beta_V3)|sumc(grad_g1_V3)|sumc(grad_g2_V3)|sumc(grad_g3_V3));
                endif;
                
            elseif capJ eq 5;
                
                if kv eq 1;
                    score_anal 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_g4)|sumc(grad_beta_V1)|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1)|sumc(grad_g4_V1));
                elseif kv eq 2;
                    score_anal 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_g4)|sumc(grad_beta_V1)|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1)|sumc(grad_g4_V1)|sumc(grad_beta_V2)|sumc(grad_g1_V2)|sumc(grad_g2_V2)|sumc(grad_g3_V2)|sumc(grad_g4_V2));
                elseif kv eq 3;
                    score_anal 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_g4)|sumc(grad_beta_V1)|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1)|sumc(grad_g4_V1)|sumc(grad_beta_V2)|sumc(grad_g1_V2)|sumc(grad_g2_V2)|sumc(grad_g3_V2)|sumc(grad_g4_V2)|
                        sumc(grad_beta_V3)|sumc(grad_g1_V3)|sumc(grad_g2_V3)|sumc(grad_g3_V3)|sumc(grad_g4_V3));
                endif;
                
            endif;
            
            
            if capJ eq 3;
                
                if kv eq 1;
                    score_anal_VE 	= (sumc(grad_beta)|sumc(grad_g1+grad_g1_V1)|sumc(grad_g2+grad_g2_V1)|sumc(grad_beta_V1));
                elseif kv eq 2;
                    score_anal_VE 	= (sumc(grad_beta)|sumc(grad_g1+grad_g1_V1+grad_g1_V2)|sumc(grad_g2+grad_g2_V1+grad_g2_V2)|sumc(grad_beta_V1)|sumc(grad_beta_V2));
                elseif kv eq 3;
                    score_anal_VE 	= (sumc(grad_beta)|sumc(grad_g1+grad_g1_V1+grad_g1_V2+grad_g1_V3)|sumc(grad_g2+grad_g2_V1+grad_g2_V2+grad_g2_V3)|
                        sumc(grad_beta_V1)|sumc(grad_beta_V2)|sumc(grad_beta_V3));
                endif;
                
            elseif capJ eq 4;
                
                if kv eq 1;
                    score_anal_VE 	= (sumc(grad_beta)|sumc(grad_g1+grad_g1_V1)|sumc(grad_g2+grad_g2_V1)|sumc(grad_g3+grad_g3_V1)|sumc(grad_beta_V1));
                elseif kv eq 2;
                    score_anal_VE 	= (sumc(grad_beta)|sumc(grad_g1+grad_g1_V1+grad_g1_V2)|sumc(grad_g2+grad_g2_V1+grad_g2_V2)|sumc(grad_g3+grad_g3_V1+grad_g3_V2)|sumc(grad_beta_V1)|sumc(grad_beta_V2));
                elseif kv eq 3;
                    score_anal_VE 	= (sumc(grad_beta)|sumc(grad_g1+grad_g1_V1+grad_g1_V2+grad_g1_V3)|sumc(grad_g2+grad_g2_V1+grad_g2_V2+grad_g2_V3)|sumc(grad_g3+grad_g3_V1+grad_g3_V2+grad_g3_V3)|
                        sumc(grad_beta_V1)|sumc(grad_beta_V2)|sumc(grad_beta_V3));
                endif;
                
            elseif capJ eq 5;
                
                if kv eq 1;
                    score_anal_VE 	= (sumc(grad_beta)|sumc(grad_g1+grad_g1_V1)|sumc(grad_g2+grad_g2_V1)|sumc(grad_g3+grad_g3_V1)|sumc(grad_g4+grad_g4_V1)|sumc(grad_beta_V1));
                elseif kv eq 2;
                    score_anal_VE 	= (sumc(grad_beta)|sumc(grad_g1+grad_g1_V1+grad_g1_V2)|sumc(grad_g2+grad_g2_V1+grad_g2_V2)|sumc(grad_g3+grad_g3_V1+grad_g3_V2)|sumc(grad_g4+grad_g4_V1+grad_g4_V2)|sumc(grad_beta_V1)|sumc(grad_beta_V2));
                elseif kv eq 3;
                    score_anal_VE 	= (sumc(grad_beta)|sumc(grad_g1+grad_g1_V1+grad_g1_V2+grad_g1_V3)|sumc(grad_g2+grad_g2_V1+grad_g2_V2+grad_g2_V3)|sumc(grad_g3+grad_g3_V1+grad_g3_V2+grad_g3_V3)|sumc(grad_g4+grad_g4_V1+grad_g4_V2+grad_g4_V3)|
                        sumc(grad_beta_V1)|sumc(grad_beta_V2)|sumc(grad_beta_V3));
                endif;
                
            endif;
            
            if capJ eq 3;
                
                if kv eq 1;
                    score_anal_RC 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_beta_V1[.,1])|sumc(grad_g1_V1)|sumc(grad_g2_V1));
                elseif kv eq 2;
                    score_anal_RC 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_beta_V1[.,1])|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_beta_V2[.,1])|sumc(grad_g1_V2)|sumc(grad_g2_V2));
                elseif kv eq 3;
                    score_anal_RC 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_beta_V1[.,1])|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_beta_V2[.,1])|sumc(grad_g1_V2)|sumc(grad_g2_V2)|
                        sumc(grad_beta_V3[.,1])|sumc(grad_g1_V3)|sumc(grad_g2_V3));
                endif;
                
            elseif capJ eq 4;
                
                if kv eq 1;
                    score_anal_RC 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_beta_V1[.,1])|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1));
                elseif kv eq 2;
                    score_anal_RC 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_beta_V1[.,1])|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1)|sumc(grad_beta_V2[.,1])|sumc(grad_g1_V2)|sumc(grad_g2_V2)|sumc(grad_g3_V2));
                elseif kv eq 3;
                    score_anal_RC 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_beta_V1[.,1])|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1)|sumc(grad_beta_V2[.,1])|sumc(grad_g1_V2)|sumc(grad_g2_V2)|sumc(grad_g3_V2)|
                        sumc(grad_beta_V3[.,1])|sumc(grad_g1_V3)|sumc(grad_g2_V3)|sumc(grad_g3_V3));
                endif;
                
            elseif capJ eq 5;
                
                if kv eq 1;
                    score_anal_RC 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_g4)|sumc(grad_beta_V1[.,1])|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1)|sumc(grad_g4_V1));
                elseif kv eq 2;
                    score_anal_RC 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_g4)|sumc(grad_beta_V1[.,1])|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1)|sumc(grad_g4_V1)|sumc(grad_beta_V2[.,1])|sumc(grad_g1_V2)|sumc(grad_g2_V2)|sumc(grad_g3_V2)|sumc(grad_g4_V2));
                elseif kv eq 3;
                    score_anal_RC 	= (sumc(grad_beta)|sumc(grad_g1)|sumc(grad_g2)|sumc(grad_g3)|sumc(grad_g4)|sumc(grad_beta_V1[.,1])|sumc(grad_g1_V1)|sumc(grad_g2_V1)|sumc(grad_g3_V1)|sumc(grad_g4_V1)|sumc(grad_beta_V2[.,1])|sumc(grad_g1_V2)|sumc(grad_g2_V2)|sumc(grad_g3_V2)|sumc(grad_g4_V2)|
                        sumc(grad_beta_V3[.,1])|sumc(grad_g1_V3)|sumc(grad_g2_V3)|sumc(grad_g3_V3)|sumc(grad_g4_V3));
                endif;
                
            endif;
            
            if sig_restrict eq 0;
                
                Score_anal      = Score_anal|sumc(grad_theta);
                Score_anal_RC   = Score_anal_RC|sumc(grad_theta);
                Score_anal_VE   = Score_anal_VE|sumc(grad_theta);
                
            endif;
            
            trap 1;
            
            if ismiss(eye(rows(bhhh_anal)) / bhhh_anal);
                LM_test     = (-999~-999~-999);
            else;
                
                V_inv     = eye(rows(bhhh_anal)) / bhhh_anal;
                LM_test 	= score_anal'*(V_inv)*score_anal;
                LM_test		= LM_test~df~cdfChic(LM_test,df);
                if LM_test[1] le 0;
                    LM_test = -999~-999~-999;
                endif;
                
                if print_it eq 1;
                    
                    ?;
                    "average score        ";;
                    meanc(abs(sumc(score_anal)));
                    "LM_test                ";;
                    LM_test;
                    ?;
                    
                endif;
                
                if check_con eq 1 and LM_test[1,3] lt 0.05;
                    
                    if count_hop lt 1;
                        
                        if sig_restrict eq 1;
                            start_HOP = zeros(kx_struct+kv+kx_mu-1+(kx_mu*(capJ-2)),1);
                        else;
                            start_HOP = zeros(kx_struct+kv+kx_mu-1+(kx_mu*(capJ-2)),1)|sig_V_inv;
                        endif;
                        count_hop = count_hop + 1;
                        goto next_try;
                        
                    endif;
                    
                endif;
                
            endif;
            
            if ismiss(eye(rows(bhhh_anal_VE)) / bhhh_anal_VE);
                LM_test_VE 	= -999~-999~-999;
            else;
                
                V_inv_VE     = eye(rows(bhhh_anal_VE)) / bhhh_anal_VE;
                LM_test_VE 	= score_anal_VE'*(V_inv_VE)*score_anal_VE;
                LM_test_VE	= LM_test_VE~df_VE~cdfChic(LM_test_VE,df_VE);
                if LM_test_VE[1] le 0;
                    LM_test_VE = -999~-999~-999;
                endif;
                
                if print_it eq 1;
                    
                    ?;
                    "average score VE         ";;
                    meanc(abs(score_anal_VE));
                    "LM_test_VE                 ";;
                    LM_test_VE;
                    ?;
                    
                endif;
                
            endif;
            
            
            if ismiss(eye(rows(bhhh_anal_RC)) / bhhh_anal_RC);
                LM_test_RC  = -999~-999~-999;
                
            else;
                
                V_inv_RC     = eye(rows(bhhh_anal_RC)) / bhhh_anal_RC;
                LM_test_RC 	= score_anal_RC'*(V_inv_RC)*score_anal_RC;
                LM_test_RC	= LM_test_RC~df_RC~cdfChic(LM_test_RC,df_RC);
                if LM_test_RC[1] le 0;
                    LM_test_RC = -999~-999~-999;
                endif;
                
                if print_it eq 1;
                    
                    ?;
                    "average score RC         ";;
                    meanc(abs(score_anal_RC));
                    "LM_test_RC                 ";;
                    LM_test_RC;
                    ?;
                    
                endif;
                
            endif;
            
            if cond_Y_N eq 1;
                
                if abs(LM_test[1,3]) le 0.05;
                    
                    if cond(chol_bhhh) ge 1500;
                        
                        condition   = eye(rows(bhhh_anal)) .* 0.00001;
                        
                        if LM_test[1,3] le 0.01;
                            condition   = eye(rows(bhhh_anal)) .* 0.0001;
                        endif;
                        
                        cond_bhhh   = bhhh_anal + condition;
                        chol_bhhh	= chol(cond_bhhh);
                        V_inv		= (eye(rows(chol_bhhh))/chol_bhhh)*((eye(rows(chol_bhhh))/chol_bhhh))';
                        LM_cond     = score_anal'*(V_inv ./ N)*score_anal;
                        LM_cond		= LM_cond~df~cdfChic(LM_cond,df);
                        
                        if LM_test[1] ge LM_cond[1];
                            LM_test = LM_cond;
                        endif;
                        
                        "LM_test (post conditioning)";;
                        LM_test;
                        ?;
                        
                    endif;
                    
                endif;
                
                if abs(LM_test_VE[1,3]) le 0.05;
                    
                    if cond(chol_bhhh_VE) ge 1200;
                        
                        condition   = eye(rows(bhhh_anal_VE)) .* 0.00001;
                        
                        if LM_test_VE[1,3] le 0.01;
                            condition   = eye(rows(bhhh_anal_VE)) .* 0.0001;
                        endif;
                        
                        cond_bhhh_VE    = bhhh_anal_VE + condition;
                        chol_bhhh_VE	= chol(cond_bhhh_VE);
                        V_inv_VE		= (eye(rows(chol_bhhh_VE))/chol_bhhh_VE)*((eye(rows(chol_bhhh_VE))/chol_bhhh_VE))';
                        LM_cond_VE      = score_anal_VE'*(V_inv_VE ./ N)*score_anal_VE;
                        LM_cond_VE		= LM_cond_VE~df_VE~cdfChic(LM_cond_VE,df_VE);
                        
                        if LM_test_VE[1] ge LM_cond_VE[1];
                            LM_test_VE = LM_cond_VE;
                        endif;
                        
                        "LM_test VE (post conditioning)";;
                        LM_test_VE;
                        ?;
                        
                    endif;
                    
                endif;
                
                if abs(LM_test_RC[1,3]) le 0.05;
                    
                    if cond(chol_bhhh_RC) ge 5000;
                        
                        condition   = eye(rows(bhhh_anal_RC)) .* 0.00001;
                        
                        if LM_test_RC[1,3] le 0.01;
                            condition   = eye(rows(bhhh_anal_RC)) .* 0.0001;
                        endif;
                        
                        cond_bhhh_RC    = bhhh_anal_RC + condition;
                        chol_bhhh_RC	= chol(cond_bhhh_RC);
                        V_inv_RC		= (eye(rows(chol_bhhh_RC))/chol_bhhh_RC)*((eye(rows(chol_bhhh_RC))/chol_bhhh_RC))';
                        LM_cond_RC      = score_anal_RC'*(V_inv_RC ./ N)*score_anal_RC;
                        LM_cond_RC		= LM_cond_RC~df_RC~cdfChic(LM_cond_RC,df_RC);
                        
                        if LM_test_RC[1] ge LM_cond_RC[1];
                            LM_test_RC = LM_cond_RC;
                        endif;
                        
                        "LM_test RC (post conditioning)";;
                        LM_test_RC;
                        ?;
                        
                    endif;
                    
                endif;
                
            endif;
            
            LM_test = LM_test|LM_test_VE|LM_test_RC;
            
            skip_LM:
                
                if LR_Yes_No eq 2;
                
                Wald_test   = -999~-999~-999;
                LR_test     = Wald_test;
                goto skip_test;
                
            endif;
            
            skip_test:
                
            endif;
        
    endif;
    
    if mtype eq 2;
        goto finish;
    endif;
    
    finish:
        
        /* NEED TO STACK IC'S TO DETERMINE OPTIMAL MODEL: KEEP HOPIT VARIANTS SEPARATE */
        
        if mtype eq 1;
        IC_ALL		= IC_OP;
    endif;
    
    if mtype eq 2;
        IC_HOPITS	= IC_HOPIT;
        IC_ALL		= IC_OP;
    endif;
    
    IC_ALL		= seqa(1,1,rows(IC_ALL))~IC_ALL;
    IC_HOPITS	= seqa(1,1,rows(IC_HOPITS))~IC_HOPITS;
    
    /* RETURN OUTPUTS */
    
    retp(OLS_B,
        OP_B, me_OP, me_OPse, V_OP, L_OP, IC_OP, pordp, y_star_OP,
        HOPIT_B, me_HOPIT, me_HOPITse, V_HOPIT, L_HOPIT, IC_HOPIT, p_HOPIT, y_star_HOPIT,
        IC_ALL, IC_HOPITS, LM_test);
endp;


/* ORDERED PROBIT MODEL */

proc ORDERED_MLE(struct PV p, struct DS d, ind);
    local x, k, d_ij, beta_Y, mu, xbeta, p_ij;
    local mu_mat, mu_xbeta;
    local grad_beta, grad_betaJ, grad_mu, grad_muJ, grad_mu1, grad_mu2, grad_mu3, grad_mu4, L_i;
    
    
    struct modelResults mm;
    
    /* DEFINE PARAMETERS */
    
    beta_Y  = pvUnpack(p,"beta_Y");
    mu      = pvUnpack(p,"mu");
    
    
    /* DEFINE DATA */
    
    d_ij	= d[1].dataMatrix;
    x 		= d[2].dataMatrix;
    
    
    /* DETERMINISTIC PART & PROBABILITIES */
    
    xbeta 	= x*beta_Y;
    mu_mat	= mu' .* ones(rows(d_ij),capJ-1);
    
    mu_xbeta	= mu_mat - (xbeta .*. ones(1,capJ-1));
    p_ij		= cdfn(mu_xbeta);
    p_ij		= p_ij[.,1]~(p_ij[.,2:cols(p_ij)] - p_ij[.,1:cols(p_ij)-1]);
    
    p_ij		= p_ij~(1-sumr(p_ij));
    
    L_i         = (sumr(d_ij .* p_ij));
    L_i         = L_i + ((L_i .lt (1e-8)) .* (1e-8));
    
    mm.Function = ln(L_i);
    
    if ind[2];
        
        /* beta_Y derivatives */
        
        {grad_beta} = _grad_beta_OP(mu_xbeta,d_ij,p_ij,x);
        
        /* mu derivatives */
        
        {grad_mu} = _grad_mu_OP(mu_xbeta,d_ij,p_ij);
        
        mm.Gradient 	= grad_beta~grad_mu;
        
    endif;
    
    
    retp(mm);
endp;

/* HOPIT/VIGNETTES MODEL */

proc HOPIT_MLE(struct PV p, struct DS d, ind);
    local X_struct, X_mu, V_x, V_ij1, V_ij2, V_ij3, V_ij4, V_ij5, jrep, d_ij, betaX, V_cons, sig_V_inv, gama, gama0, xbeta, zgama, p_ij;
    local mu_ij, mu_ijV, mu_mat, mu_xbeta, lnLiHOPIT;
    local mu_xbetaV1, p_ijV1, lnLiV1, mu_xbetaV2, p_ijV2, lnLiV2, mu_xbetaV3, p_ijV3, lnLiV3, mu_xbetaV4, p_ijV4, lnLiV4, mu_xbetaV5, p_ijV5, lnLiV5;
    local LnLiV, sum_check, sum_check1, sum_check2, sum_check3;
    local grad_beta, grad_V_con, grad_V_con1, grad_V_con2, grad_V_con3, grad_g1V1, grad_g1V2, grad_g1V3, grad_g2V1, grad_g2V2, grad_g2V3;
    local grad_theta, grad_theta_V1, grad_theta_V2, grad_theta_V3, grad_g3V1, grad_g3V2, grad_g3V3, grad_g4V1, grad_g4V2, grad_g4V3;
    local grad_g1, grad_g2, grad_g3, grad_g4, grad_betaJ;
    local grad_g1a, grad_g1b, grad_g1V1a, grad_g1V1b, grad_g1V2a, grad_g1V2b, grad_g1V3a, grad_g1V3b;
    local grad_g2a, grad_g2b, grad_g2V1a, grad_g2V1b, grad_g2V2a, grad_g2V2b, grad_g2V3a, grad_g2V3b;
    local grad_g3a, grad_g3b, grad_g3V1a, grad_g3V1b, grad_g3V2a, grad_g3V2b, grad_g3V3a, grad_g3V3b;
    local grad_g4a, grad_g4b, grad_g4V1a, grad_g4V1b, grad_g4V2a, grad_g4V2b, grad_g4V3a, grad_g4V3b;
    
    clear lnLiV1, lnLiV2, lnLiV3, lnLiV4, lnLiV5;
    
    struct modelResults mm;
    
    /* DEFINE PARAMETERS */
    
    betaX    	= pvUnpack(p,"beta_Y");
    V_cons		= pvUnpack(p,"V_cons");
    gama		= pvUnpack(p,"gama");
    gama0		= pvUnpack(p,"gama0")';
    gama        = (0~gama0)|gama;
    gama	    = gama';
    if sig_restrict eq 0;
        sig_V_inv	= pvUnpack(p,"sig_V_inv");
    else;
        sig_V_inv	= 1;
    endif;
    
    /* DEFINE DATA */
    
    d_ij		= d[1].dataMatrix;
    X_struct 	= d[2].dataMatrix;
    X_mu		= d[3].dataMatrix;
    
    
    V_x			= X_mu;
    V_ij1		= d[4].dataMatrix;
    
    if kv ge 2;
        V_ij2		= d[5].dataMatrix;
        if kv ge 3;
            V_ij3	= d[6].dataMatrix;
        endif;
    endif;
    
    /* DEFINE COMMON BOUNDARY PARAMETERS */
    
    if exp_mu eq 1;
        
        mu_ij	= _mu_ij(X_mu*gama);
        mu_ijV	= _mu_ij(V_x*gama);
        
    elseif exp_mu eq 0;
        
        mu_ij	= X_mu*gama;
        mu_ijV	= V_x*gama;
        
    endif;
    
    /* HOPIT FOR THE VIGNETTES */
    
    mu_xbetaV1	= (mu_ijV - V_cons[1]) .* sig_V_inv;
    p_ijV1		= cdfn(mu_xbetaV1);
    p_ijV1		= p_ijV1[.,1]~(p_ijV1[.,2:cols(p_ijV1)] - p_ijV1[.,1:cols(p_ijV1)-1]);
    p_ijV1		= p_ijV1~(1-sumr(p_ijV1));
    //lnLiV1		= ln(sumr(V_ij1 .* p_ijV1));
    
    lnLiV1   = (sumr(V_ij1 .* p_ijV1));
    //lnLiV1   = lnLiV1 + ((lnLiV1 .lt (1e-8)) .* (1e-8));
    lnLiV1   = ln(lnLiV1);
    
    if maxc(maxc(p_ijV1)) .gt 1;
        
        "P>=1 V1 stop";
        maxc(p_ijV1)';
        ?;
        sig_V_inv;
        ?;
        meanc(mu_ij)';
        stop;
        
    endif;
    
    if flt(minc(vec(p_ijV1)), 0);
        "minc(vec(p_ijV1))";;
        minc(vec(p_ijV1));
        stop;
    endif;
    
    if kv ge 2;
        
        mu_xbetaV2	= (mu_ijV - V_cons[2]) .* sig_V_inv;
        p_ijV2		= cdfn(mu_xbetaV2);
        p_ijV2		= p_ijV2[.,1]~(p_ijV2[.,2:cols(p_ijV2)] - p_ijV2[.,1:cols(p_ijV2)-1]);
        p_ijV2		= p_ijV2~(1-sumr(p_ijV2));
        //lnLiV2		= ln(sumr(V_ij2 .* p_ijV2));
        
        lnLiV2   = (sumr(V_ij2 .* p_ijV2));
        //lnLiV2   = lnLiV2 + ((lnLiV2 .lt (1e-8)) .* (1e-8));
        lnLiV2   = ln(lnLiV2);
        
        if maxc(maxc(p_ijV2)) .gt 1;
            "P>=1 V2 stop";
            maxc(p_ijV2)';
            stop;
        endif;
        
        if flt(minc(vec(p_ijV2)), 0);
            "minc(vec(p_ijV2))";;
            minc(vec(p_ijV2));
            stop;
        endif;
        
        if kv ge 3;
            
            mu_xbetaV3	= (mu_ijV - V_cons[3]) .* sig_V_inv;
            p_ijV3		= cdfn(mu_xbetaV3);
            p_ijV3		= p_ijV3[.,1]~(p_ijV3[.,2:cols(p_ijV3)] - p_ijV3[.,1:cols(p_ijV3)-1]);
            p_ijV3		= p_ijV3~(1-sumr(p_ijV3));
            //lnLiV3		= ln(sumr(V_ij3 .* p_ijV3));
            
            lnLiV3   = (sumr(V_ij3 .* p_ijV3));
            lnLiV3   = lnLiV3 + ((lnLiV3 .lt (1e-8)) .* (1e-8));
            lnLiV3   = ln(lnLiV3);
            
            if maxc(maxc(p_ijV3)) .gt 1;
                "P>=1 V3 stop";
                maxc(p_ijV3)';
                stop;
            endif;
            
            if flt(minc(vec(p_ijV3)), 0);
                "minc(vec(p_ijV3))";;
                minc(vec(p_ijV3));
                stop;
            endif;
            
        endif;
        
    endif;
    
    if kv eq 1;
        LnLiV = (lnLiV1);
    elseif kv eq 2;
        LnLiV = (lnLiV1) + (lnLiV2);
    elseif kv eq 3;
        LnLiV = (lnLiV1) + (lnLiV2) + (lnLiV3);
    endif;
    
    /* DETERMINISTIC PART & PROBABILITIES */
    
    xbeta 		= X_struct*betaX;
    mu_xbeta	= mu_ij - (xbeta .*. ones(1,capJ-1));
    p_ij		= cdfn(mu_xbeta);
    p_ij		= p_ij[.,1]~(p_ij[.,2:cols(p_ij)] - p_ij[.,1:cols(p_ij)-1]);
    p_ij		= p_ij~(1-sumr(p_ij));
    
    lnLiHOPIT   = (sumr(d_ij .* p_ij));
    lnLiHOPIT   = lnLiHOPIT + ((lnLiHOPIT .lt (1e-8)) .* (1e-8));
    lnLiHOPIT   = ln(lnLiHOPIT);
    
    if flt(minc(vec(p_ij)), 0);
        "minc(vec(p_ij))";;
        minc(vec(p_ij));
        stop;
    endif;
    
    if ind[1];
        mm.Function = lnLiV + lnLiHOPIT;
    endif;
    
    
    if ind[2];
        
        /* beta_Y derivatives */
        
        {grad_beta} = _grad_beta_OP(mu_xbeta,d_ij,p_ij,X_struct);
        
        
        /*  vignette constants derivatives */
        
        if kv eq 1;
            
            {grad_V_con1} 	= _grad_V_HOPIT(sig_V_inv,mu_xbetaV1,p_ijV1,V_ij1);
            grad_V_con      = grad_V_con1;
            
        elseif kv eq 2;
            
            {grad_V_con1} 	= _grad_V_HOPIT(sig_V_inv,mu_xbetaV1,p_ijV1,V_ij1);
            {grad_V_con2} 	= _grad_V_HOPIT(sig_V_inv,mu_xbetaV2,p_ijV2,V_ij2);
            grad_V_con      = grad_V_con1~grad_V_con2;
            
        elseif kv eq 3;
            
            {grad_V_con1} 	= _grad_V_HOPIT(sig_V_inv,mu_xbetaV1,p_ijV1,V_ij1);
            {grad_V_con2} 	= _grad_V_HOPIT(sig_V_inv,mu_xbetaV2,p_ijV2,V_ij2);
            {grad_V_con3} 	= _grad_V_HOPIT(sig_V_inv,mu_xbetaV3,p_ijV3,V_ij3);
            grad_V_con      = grad_V_con1~grad_V_con2~grad_V_con3;
            
        endif;
        
        if exp_mu eq 1;
            
            /* gradient for gamma_1 HOPIT */
            
            {grad_g1}		= _grad_g_j(gama[2:kx_mu,1],X_mu[.,2:cols(X_mu)],mu_xbeta,p_ij,d_ij,1);
            
            /* gradient for gamma_1 V1 */
            
            {grad_g1V1}		= _grad_g_j(gama[2:kx_mu,1],V_x[.,2:cols(V_x)],mu_xbetaV1,p_ijV1,V_ij1,sig_V_inv);
            
            if kv ge 2;
                
                /* gradient for gamma_1 V2 */
                
                {grad_g1V2}		= _grad_g_j(gama[2:kx_mu,1],V_x[.,2:cols(V_x)],mu_xbetaV2,p_ijV2,V_ij2,sig_V_inv);
                
            endif;
            
            if kv ge 3;
                
                /* gradient for gamma_1 V3 */
                
                {grad_g1V3}		= _grad_g_j(gama[2:kx_mu,1],V_x[.,2:cols(V_x)],mu_xbetaV3,p_ijV3,V_ij3,sig_V_inv);
                
            endif;
            
            if kv eq 1;
                grad_g1 	= grad_g1 + grad_g1V1;
            elseif kv eq 2;
                grad_g1 	= grad_g1 + grad_g1V1 + grad_g1V2;
            elseif kv eq 3;
                grad_g1 	= grad_g1 + grad_g1V1 + grad_g1V2 + grad_g1V3;
            endif;
            
            if capJ gt 3;
                
                /* gradient for gamma_2 HOPIT */
                
                {grad_g2}		= _grad_g_j(gama[.,2],X_mu,mu_xbeta[.,2:cols(mu_xbeta)],p_ij[.,2:capJ],d_ij[.,2:capJ],1);
                
                /* gradient for gamma_2 V1 */
                
                {grad_g2V1}		= _grad_g_j(gama[.,2],V_x,mu_xbetaV1[.,2:cols(mu_xbetaV1)],p_ijV1[.,2:capJ],V_ij1[.,2:capJ],sig_V_inv);
                
                if kv ge 2;
                    
                    /* gradient for gamma_2 V2 */
                    
                    {grad_g2V2}		= _grad_g_j(gama[.,2],V_x,mu_xbetaV2[.,2:cols(mu_xbetaV2)],p_ijV2[.,2:capJ],V_ij2[.,2:capJ],sig_V_inv);
                    
                endif;
                
                if kv ge 3;
                    
                    /* gradient for gamma_2 V3 */
                    
                    {grad_g2V3}		= _grad_g_j(gama[.,2],V_x,mu_xbetaV3[.,2:cols(mu_xbetaV3)],p_ijV3[.,2:capJ],V_ij3[.,2:capJ],sig_V_inv);
                    
                endif;
                
            elseif capJ eq 3;
                
                {grad_g2}		= _grad_gJ(gama[.,2],X_mu,mu_xbeta[.,cols(mu_xbeta)],p_ij[.,2:capJ],d_ij[.,2:capJ],1);
                
                /* gradient for gamma_2 V1 */
                
                {grad_g2V1}		= _grad_gJ(gama[.,2],V_x,mu_xbetaV1[.,cols(mu_xbetaV1)],p_ijV1[.,2:capJ],V_ij1[.,2:capJ],sig_V_inv);
                
                if kv ge 2;
                    
                    /* gradient for gamma_2 V2 */
                    
                    {grad_g2V2}		= _grad_gJ(gama[.,2],V_x,mu_xbetaV2[.,cols(mu_xbetaV2)],p_ijV2[.,2:capJ],V_ij2[.,2:capJ],sig_V_inv);
                    
                endif;
                
                if kv ge 3;
                    
                    /* gradient for gamma_2 V3 */
                    
                    {grad_g2V3}		= _grad_gJ(gama[.,2],V_x,mu_xbetaV3[.,cols(mu_xbetaV3)],p_ijV3[.,2:capJ],V_ij3[.,2:capJ],sig_V_inv);
                    
                endif;
                
            endif;
            
            if kv eq 1;
                grad_g2 	= grad_g2 + grad_g2V1;
            elseif kv eq 2;
                grad_g2 	= grad_g2 + grad_g2V1 + grad_g2V2;
            elseif kv eq 3;
                grad_g2 	= grad_g2 + grad_g2V1 + grad_g2V2 + grad_g2V3;
            endif;
            
            if capJ ge 5;
                
                /* gradient for gamma_3 HOPIT */
                
                {grad_g3}		= _grad_g_j(gama[.,3],X_mu,mu_xbeta[.,3:cols(mu_xbeta)],p_ij[.,3:capJ],d_ij[.,3:capJ],1);
                
                /* gradient for gamma_3 V1 */
                
                {grad_g3V1}		= _grad_g_j(gama[.,3],V_x,mu_xbetaV1[.,3:cols(mu_xbetaV1)],p_ijV1[.,3:capJ],V_ij1[.,3:capJ],sig_V_inv);
                
                if kv ge 2;
                    
                    /* gradient for gamma_3 V2 */
                    
                    {grad_g3V2}		= _grad_g_j(gama[.,3],V_x,mu_xbetaV2[.,3:cols(mu_xbetaV2)],p_ijV2[.,3:capJ],V_ij2[.,3:capJ],sig_V_inv);
                    
                endif;
                
                if kv ge 3;
                    
                    /* gradient for gamma_3 V3 */
                    
                    {grad_g3V3}		= _grad_g_j(gama[.,3],V_x,mu_xbetaV3[.,3:cols(mu_xbetaV3)],p_ijV3[.,3:capJ],V_ij3[.,3:capJ],sig_V_inv);
                    
                endif;
                
            elseif capJ eq 4;
                
                {grad_g3}		= _grad_gJ(gama[.,3],X_mu,mu_xbeta[.,cols(mu_xbeta)],p_ij[.,3:capJ],d_ij[.,3:capJ],1);
                
                /* gradient for gamma_3 V1 */
                
                {grad_g3V1}		= _grad_gJ(gama[.,3],V_x,mu_xbetaV1[.,cols(mu_xbetaV1)],p_ijV1[.,3:capJ],V_ij1[.,3:capJ],sig_V_inv);
                
                if kv ge 2;
                    
                    /* gradient for gamma_3 V2 */
                    
                    {grad_g3V2}		= _grad_gJ(gama[.,3],V_x,mu_xbetaV2[.,cols(mu_xbetaV2)],p_ijV2[.,3:capJ],V_ij2[.,3:capJ],sig_V_inv);
                    
                endif;
                
                if kv ge 3;
                    
                    /* gradient for gamma_3 V3 */
                    
                    {grad_g3V3}		= _grad_gJ(gama[.,3],V_x,mu_xbetaV3[.,cols(mu_xbetaV3)],p_ijV3[.,3:capJ],V_ij3[.,3:capJ],sig_V_inv);
                    
                endif;
                
            endif;
            
            if capJ gt 3;
                
                if kv eq 1;
                    grad_g3 		= grad_g3 + grad_g3V1;
                elseif kv eq 2;
                    grad_g3 		= grad_g3 + grad_g3V1 + grad_g3V2;
                elseif kv eq 3;
                    grad_g3 		= grad_g3 + grad_g3V1 + grad_g3V2 + grad_g3V3;
                endif;
                
            endif;
            
            if capJ ge 5;
                
                /* gradient for gamma_4 HOPIT */
                
                {grad_g4}		= _grad_gJ(gama[.,4],X_mu,mu_xbeta[.,cols(mu_xbeta)],p_ij[.,4:capJ],d_ij[.,4:capJ],1);
                
                /* gradient for gamma_4 V1 */
                
                {grad_g4V1}		= _grad_gJ(gama[.,4],V_x,mu_xbetaV1[.,cols(mu_xbetaV1)],p_ijV1[.,4:capJ],V_ij1[.,4:capJ],sig_V_inv);
                
                if kv ge 2;
                    
                    /* gradient for gamma_4 V2 */
                    
                    {grad_g4V2}		= _grad_gJ(gama[.,4],V_x,mu_xbetaV2[.,cols(mu_xbetaV2)],p_ijV2[.,4:capJ],V_ij2[.,4:capJ],sig_V_inv);
                    
                endif;
                
                if kv ge 3;
                    
                    /* gradient for gamma_4 V3 */
                    
                    {grad_g4V3}		= _grad_gJ(gama[.,4],V_x,mu_xbetaV3[.,cols(mu_xbetaV3)],p_ijV3[.,4:capJ],V_ij3[.,4:capJ],sig_V_inv);
                    
                endif;
                
                if kv eq 1;
                    grad_g4 		= grad_g4 + grad_g4V1;
                elseif kv eq 2;
                    grad_g4 		= grad_g4 + grad_g4V1 + grad_g4V2;
                elseif kv eq 3;
                    grad_g4 		= grad_g4 + grad_g4V1 + grad_g4V2 + grad_g4V3;
                endif;
                
            endif;
            
        elseif exp_mu eq 0;
            
            /* gradient for gamma_1 HOPIT */
            
            {grad_g1a}		= _grad_g_j_linear(X_mu[.,2:kx_mu],mu_xbeta[.,1],p_ij[.,1],d_ij[.,1],1);
            {grad_g1b}		= _grad_g_j_linear(X_mu[.,2:kx_mu],mu_xbeta[.,1],p_ij[.,2],d_ij[.,2],1);
            grad_g1			= grad_g1a - grad_g1b;
            
            /* gradient for gamma_1 V1 */
            
            {grad_g1V1a}		= _grad_g_j_linear(V_x[.,2:kx_mu],mu_xbetaV1[.,1],p_ijV1[.,1],V_ij1[.,1],sig_V_inv);
            {grad_g1V1b}		= _grad_g_j_linear(V_x[.,2:kx_mu],mu_xbetaV1[.,1],p_ijV1[.,2],V_ij1[.,2],sig_V_inv);
            grad_g1V1			= grad_g1V1a - grad_g1V1b;
            
            if kv ge 2;
                
                /* gradient for gamma_1 V2 */
                
                {grad_g1V2a}		= _grad_g_j_linear(V_x[.,2:kx_mu],mu_xbetaV2[.,1],p_ijV2[.,1],V_ij2[.,1],sig_V_inv);
                {grad_g1V2b}		= _grad_g_j_linear(V_x[.,2:kx_mu],mu_xbetaV2[.,1],p_ijV2[.,2],V_ij2[.,2],sig_V_inv);
                grad_g1V2			= grad_g1V2a - grad_g1V2b;
                
            endif;
            
            if kv ge 3;
                
                /* gradient for gamma_1 V3 */
                
                {grad_g1V3a}		= _grad_g_j_linear(V_x[.,2:kx_mu],mu_xbetaV3[.,1],p_ijV3[.,1],V_ij3[.,1],sig_V_inv);
                {grad_g1V3b}		= _grad_g_j_linear(V_x[.,2:kx_mu],mu_xbetaV3[.,1],p_ijV3[.,2],V_ij3[.,2],sig_V_inv);
                grad_g1V3			= grad_g1V3a - grad_g1V3b;
                
            endif;
            
            if kv eq 1;
                grad_g1 	= grad_g1 + grad_g1V1;
            elseif kv eq 2;
                grad_g1 	= grad_g1 + grad_g1V1 + grad_g1V2;
            elseif kv eq 3;
                grad_g1 	= grad_g1 + grad_g1V1 + grad_g1V2 + grad_g1V3;
            endif;
            
            /* gradient for gamma_2 HOPIT */
            
            {grad_g2a}		= _grad_g_j_linear(X_mu,mu_xbeta[.,2],p_ij[.,2],d_ij[.,2],1);
            {grad_g2b}		= _grad_g_j_linear(X_mu,mu_xbeta[.,2],p_ij[.,3],d_ij[.,3],1);
            grad_g2			= grad_g2a - grad_g2b;
            
            /* gradient for gamma_2 V1 */
            
            {grad_g2V1a}		= _grad_g_j_linear(V_x,mu_xbetaV1[.,2],p_ijV1[.,2],V_ij1[.,2],sig_V_inv);
            {grad_g2V1b}		= _grad_g_j_linear(V_x,mu_xbetaV1[.,2],p_ijV1[.,3],V_ij1[.,3],sig_V_inv);
            grad_g2V1			= grad_g2V1a - grad_g2V1b;
            
            if kv ge 2;
                
                /* gradient for gamma_2 V2 */
                
                {grad_g2V2a}		= _grad_g_j_linear(V_x,mu_xbetaV2[.,2],p_ijV2[.,2],V_ij2[.,2],sig_V_inv);
                {grad_g2V2b}		= _grad_g_j_linear(V_x,mu_xbetaV2[.,2],p_ijV2[.,3],V_ij2[.,3],sig_V_inv);
                grad_g2V2			= grad_g2V2a - grad_g2V2b;
                
            endif;
            
            if kv ge 3;
                
                /* gradient for gamma_2 V3 */
                
                {grad_g2V3a}		= _grad_g_j_linear(V_x,mu_xbetaV3[.,2],p_ijV3[.,2],V_ij3[.,2],sig_V_inv);
                {grad_g2V3b}		= _grad_g_j_linear(V_x,mu_xbetaV3[.,2],p_ijV3[.,3],V_ij3[.,3],sig_V_inv);
                grad_g2V3			= grad_g2V3a - grad_g2V3b;
                
            endif;
            
            if kv eq 1;
                grad_g2 	= grad_g2 + grad_g2V1;
            elseif kv eq 2;
                grad_g2 	= grad_g2 + grad_g2V1 + grad_g2V2;
            elseif kv eq 3;
                grad_g2 	= grad_g2 + grad_g2V1 + grad_g2V2 + grad_g2V3;
            endif;
            
            /* gradient for gamma_3 HOPIT */
            
            if capJ gt 3;
                
                {grad_g3a}		= _grad_g_j_linear(X_mu,mu_xbeta[.,3],p_ij[.,3],d_ij[.,3],1);
                {grad_g3b}		= _grad_g_j_linear(X_mu,mu_xbeta[.,3],p_ij[.,4],d_ij[.,4],1);
                grad_g3			= grad_g3a - grad_g3b;
                
                /* gradient for gamma_3 V1 */
                
                {grad_g3V1a}		= _grad_g_j_linear(V_x,mu_xbetaV1[.,3],p_ijV1[.,3],V_ij1[.,3],sig_V_inv);
                {grad_g3V1b}		= _grad_g_j_linear(V_x,mu_xbetaV1[.,3],p_ijV1[.,4],V_ij1[.,4],sig_V_inv);
                grad_g3V1			= grad_g3V1a - grad_g3V1b;
                
                if kv ge 2;
                    
                    /* gradient for gamma_3 V2 */
                    
                    {grad_g3V2a}		= _grad_g_j_linear(V_x,mu_xbetaV2[.,3],p_ijV2[.,3],V_ij2[.,3],sig_V_inv);
                    {grad_g3V2b}		= _grad_g_j_linear(V_x,mu_xbetaV2[.,3],p_ijV2[.,4],V_ij2[.,4],sig_V_inv);
                    grad_g3V2			= grad_g3V2a - grad_g3V2b;
                    
                endif;
                
                if kv ge 3;
                    
                    /* gradient for gamma_3 V3 */
                    
                    {grad_g3V3a}		= _grad_g_j_linear(V_x,mu_xbetaV3[.,3],p_ijV3[.,3],V_ij3[.,3],sig_V_inv);
                    {grad_g3V3b}		= _grad_g_j_linear(V_x,mu_xbetaV3[.,3],p_ijV3[.,4],V_ij3[.,4],sig_V_inv);
                    grad_g3V3			= grad_g3V3a - grad_g3V3b;
                    
                endif;
                
                if kv eq 1;
                    grad_g3 	= grad_g3 + grad_g3V1;
                elseif kv eq 2;
                    grad_g3 	= grad_g3 + grad_g3V1 + grad_g3V2;
                elseif kv eq 3;
                    grad_g3 	= grad_g3 + grad_g3V1 + grad_g3V2 + grad_g3V3;
                endif;
                
            endif;
            
            if capJ ge 5;
                
                /* gradient for gamma_4 HOPIT */
                
                {grad_g4a}		= _grad_g_j_linear(X_mu,mu_xbeta[.,4],p_ij[.,4],d_ij[.,4],1);
                {grad_g4b}		= _grad_g_j_linearJ(X_mu,mu_xbeta[.,4],p_ij[.,5],d_ij[.,5],1);
                grad_g4			= grad_g4a - grad_g4b;
                
                /* gradient for gamma_4 V1 */
                
                {grad_g4V1a}	= _grad_g_j_linear(V_x,mu_xbetaV1[.,4],p_ijV1[.,4],V_ij1[.,4],sig_V_inv);
                {grad_g4V1b}	= _grad_g_j_linearJ(V_x,mu_xbetaV1[.,4],p_ijV1[.,5],V_ij1[.,5],sig_V_inv);
                grad_g4V1		= grad_g4V1a - grad_g4V1b;
                
                if kv ge 2;
                    
                    /* gradient for gamma_4 V2 */
                    
                    {grad_g4V2a}	= _grad_g_j_linear(V_x,mu_xbetaV2[.,4],p_ijV2[.,4],V_ij2[.,4],sig_V_inv);
                    {grad_g4V2b}	= _grad_g_j_linearJ(V_x,mu_xbetaV2[.,4],p_ijV2[.,5],V_ij2[.,5],sig_V_inv);
                    grad_g4V2		= grad_g4V2a - grad_g4V2b;
                    
                endif;
                
                if kv eq 3;
                    
                    /* gradient for gamma_4 V3 */
                    
                    {grad_g4V3a}		= _grad_g_j_linear(V_x,mu_xbetaV3[.,4],p_ijV3[.,4],V_ij3[.,4],sig_V_inv);
                    {grad_g4V3b}		= _grad_g_j_linearJ(V_x,mu_xbetaV3[.,4],p_ijV3[.,5],V_ij3[.,5],sig_V_inv);
                    grad_g4V3			= grad_g4V3a - grad_g4V3b;
                    
                endif;
                
                if kv eq 1;
                    grad_g4 	= grad_g4 + grad_g4V1;
                elseif kv eq 2;
                    grad_g4 	= grad_g4 + grad_g4V1 + grad_g4V2;
                elseif kv eq 3;
                    grad_g4 	= grad_g4 + grad_g4V1 + grad_g4V2 + grad_g4V3;
                endif;
                
            endif;
            
        endif;
        
        /* theta = 1/sigma derivatives; V1, V2 and V3 */
        
        if sig_restrict eq 0;
            
            if kv eq 1;
                
                {grad_theta_V1}	= _grad_theta_HOPIT(mu_xbetaV1,sig_V_inv,p_ijV1,V_ij1);
                grad_theta		= grad_theta_V1;
                
            elseif kv eq 2;
                
                {grad_theta_V1}	= _grad_theta_HOPIT(mu_xbetaV1,sig_V_inv,p_ijV1,V_ij1);
                {grad_theta_V2}	= _grad_theta_HOPIT(mu_xbetaV2,sig_V_inv,p_ijV2,V_ij2);
                grad_theta		= grad_theta_V1 + grad_theta_V2;
                
                
            elseif kv eq 3;
                
                {grad_theta_V1}	= _grad_theta_HOPIT(mu_xbetaV1,sig_V_inv,p_ijV1,V_ij1);
                {grad_theta_V2}	= _grad_theta_HOPIT(mu_xbetaV2,sig_V_inv,p_ijV2,V_ij2);
                {grad_theta_V3}	= _grad_theta_HOPIT(mu_xbetaV3,sig_V_inv,p_ijV3,V_ij3);
                grad_theta		= grad_theta_V1 + grad_theta_V2 + grad_theta_V3;
                
            endif;
            
        endif;
        
        if sig_restrict eq 0;
            
            if capJ eq 5;
                mm.Gradient 	= grad_beta~grad_V_con~grad_g1~grad_g2~grad_g3~grad_g4~grad_theta;
            elseif capJ eq 4;
                mm.Gradient 	= grad_beta~grad_V_con~grad_g1~grad_g2~grad_g3~grad_theta;
            elseif capJ eq 3;
                mm.Gradient 	= grad_beta~grad_V_con~grad_g1~grad_g2~grad_theta;
            endif;
            
        else;
            
            if capJ eq 5;
                mm.Gradient 	= grad_beta~grad_V_con~grad_g1~grad_g2~grad_g3~grad_g4;
            elseif capJ eq 4;
                mm.Gradient 	= grad_beta~grad_V_con~grad_g1~grad_g2~grad_g3;
            elseif capJ eq 3;
                mm.Gradient 	= grad_beta~grad_V_con~grad_g1~grad_g2;
            endif;
            
        endif;
        
    endif;
    
    retp(mm);
endp;

/* PROBABILITY PROCEDURES FOR ME'S AND SE'S(ME'S) */

/* USE GENERIC METHOD FOR ME'S AND SE'S */

proc 4 = mes(vcov,gama,capJ,me_model,xbar);
    local G_ME, G_P, me_ses, p1, p2, h, xbargama, mes, i, g, gp;
    local prob, probability, prob_se, me0_np, me0_p, me0_se_np, me0_se_p;
    local jloop, G_MEnp, G_MEp, G_Pnp, G_Pp, mes_np, mes_p, prob_np, prob_p;
    local me_np_se, prob_np_se, me_p_se, prob_p_se, mark;
    
    clearg k, kg, mtype_me;
    
    if rows(xbar) eq 1;
        xbar = xbar';
    endif;
    
    k  = rows(xbar);
    kg = rows(gama);
    
    
    /* STACK THE XBARS AND PARAMETERS*/
    
    xbargama = xbar|gama;
    
    /* WHICH MODEL? */
    
    let mtype_me = ordered hopit;
    mtype_me = indcv(me_model,mtype_me);
    
    
    /* INITIALISE MATRICES */
    
    G_ME    = zeros(1,rows(gama));
    G_P     = zeros(1,rows(gama));
    mes     = zeros(1,1);
    prob	= zeros(1,1);
    
    G_MEnp  = zeros(1,rows(gama));
    G_MEp   = zeros(1,rows(gama));
    G_Pnp   = zeros(1,rows(gama));
    G_Pp    = zeros(1,rows(gama));
    
    clear mes,me_ses,mes_np,mes_p,me_np_se,me_p_se,prob,prob_se;
    
    
    /* CONSTRUCT P1 & P2 FOR EXTRACTING SECOND ORDER DERIVATIVES WRT PARAMETERS */
    
    p1=k+1;
    
    for i(rows(xbar)+2,rows(xbargama),1);
        p1=p1|i;
    endfor;
    
    p2=1;
    
    for i(2,rows(xbar),1);
        p2=p2|i;
    endfor;
    
    /* ESTIMATING MEs AND SEs*/
    
    if mtype_me eq 1;                              @ ORDERED PROBIT @
        
        jj          = 1;
        
        do while jj le capJ;
            
            h           = hessp(&ordp_p,xbargama);
            g           = gradp(&ordp_p,xbargama);
            probability = ordp_p(xbargama);
            gp          = gradp(&ordp_p,xbargama);
            
            h           = submat(h,p1',p2')';
            gp          = submat(gp',p1,0)';
            
            G_ME        = G_ME|h;
            G_P         = G_P|gp;
            mes         = mes|g[1:rows(xbar)]';
            prob        = prob~probability';
            
            jj          = jj+1;
            
        endo;
        
        G_ME	    = G_ME[2:rows(G_ME),.];
        G_P 	    = G_P[2:rows(G_P),.];
        me_ses	    = G_ME*vcov*G_ME';
        me_ses 	    = sqrt(diag(me_ses));
        me_ses 	    = reshape(me_ses,capJ,k)';
        mes         = mes[2:rows(mes)];
        mes         = reshape(mes,capJ,k)';
        
        prob        = prob[.,2:cols(prob)];
        prob_se	    = G_P*vcov*G_P';
        prob_se  	= sqrt(diag(prob_se));
        
    elseif mtype_me eq 2;					@ HOPIT MODEL @
        
        jj          = 1;
        
        do while jj le capJ;
            
            h           = hessp(&HOPIT_p,xbargama);
            g           = gradp(&HOPIT_p,xbargama);
            probability = HOPIT_p(xbargama);
            gp          = gradp(&HOPIT_p,xbargama);
            
            h           = submat(h,p1',p2')';
            gp          = submat(gp',p1,0)';
            
            G_ME        = G_ME|h;
            G_P         = G_P|gp;
            mes         = mes|g[1:rows(xbar)]';
            prob        = prob~probability';
            
            jj          = jj+1;
            
        endo;
        
        G_ME	    = G_ME[2:rows(G_ME),.];
        G_P 	    = G_P[2:rows(G_P),.];
        me_ses	    = G_ME*vcov*G_ME';
        me_ses 	    = sqrt(diag(me_ses));
        me_ses 	    = reshape(me_ses,capJ,k)';
        mes         = mes[2:rows(mes)];
        mes         = reshape(mes,capJ,k)';
        
        prob        = prob[.,2:cols(prob)];
        prob_se	    = G_P*vcov*G_P';
        prob_se  	= sqrt(diag(prob_se));
        
    endif;
    
    retp(mes,me_ses,
        prob,prob_se);
endp;


/* CONSTRUCT ORDERED PROBIT PROBABILITIES */

proc 1 = ordp_p(xbargama);
    local jrep, jrep2, beta_Y, mu, xbeta, mu_xbeta, p_ij, check, pord, gama;
    
    xbar        = xbargama[1:kx_struct-1];
    gama        = xbargama[kx_struct:rows(xbargama)];
    beta_Y		= gama[1:kx_struct-1];
    mu		    = gama[kx_struct:rows(gama)];
    
    if rows(gama) ne (rows(beta_Y)+rows(mu));
        "error in ordp_p(xbargama)";
        stop;
    endif;
    
    xbeta 		= xbar'*beta_Y;
    
    mu_xbeta	= mu' - (xbeta .*. ones(1,capJ-1));
    p_ij		= cdfn(mu_xbeta);
    p_ij		= p_ij[1]~(p_ij[2:cols(p_ij)] - p_ij[1:cols(p_ij)-1]);
    
    pord		= p_ij~(1-sumr(p_ij));
    
    check = meanc(sumc(pord'));
    if check gt 1.0001 or check lt .9999;
        errorlog "ERROR IN DEFINING ORDERED PROBIT PROBABILITIES; PROGRAM TERMINATED";
        "mean probablitiy =";;
        check;
        stop;
    endif;
    
    retp(pord[.,jj]);
endp;


/* CONSTRUCT HOPIT PROBABILITIES */

proc 1 = HOPIT_p(xbargama);
    local beta_Y, V_cons, sig_V, gama0, gama, xbar, params, X_struct, X_mu, zgama, mu_ij, phop, xbeta, mu_xbeta, p_ij;
    
    xbar	= xbargama[1:kx_all];
    params  = xbargama[kx_all+1:rows(xbargama)];
    beta_Y  = params[1:kx_struct];
    gama    = 0|params[kx_struct+kv+1:kx_struct+kv+kx_mu-1+(kx_mu*(capJ-2))];
    
    X_struct	= xbar'; //selif(xbar,X_struct_vec')';
    X_mu		= xbar'; //selif(xbar,X_mu_vec')';
    
    /* DEFINE COMMON BOUNDARY PARAMETERS */
    
    gama	= reshape(gama,capJ-1,kx_mu)';
    zgama	= X_mu*gama;
    if exp_mu eq 1;
        mu_ij	= _mu_ij(zgama);
    elseif exp_mu eq 0;
        mu_ij	= zgama;
    endif;
    
    /* DETERMINISTIC PART & PROBABILITIES */
    
    xbeta 		= X_struct*beta_Y;
    mu_xbeta	= mu_ij - (xbeta .*. ones(1,capJ-1));
    p_ij		= cdfn(mu_xbeta);
    p_ij		= p_ij[.,1]~(p_ij[.,2:cols(p_ij)] - p_ij[.,1:cols(p_ij)-1]);
    phop		= p_ij~(1-sumr(p_ij));
    
    retp(phop[.,jj]);
endp;


/* PROCEDURE FOR PARAMETRISING THE BOUNDRAY POINTS */

proc 1 = _mu_ij(zgama);
    local mu_ij;
    
    if exp_mu eq 1;
        
        mu_ij	= exp(zgama[.,1]);
        
        for jrep (2,capJ-1,1);
            mu_ij	= mu_ij~((mu_ij[.,jrep-1]) + exp(zgama[.,jrep]));
        endfor;
        
    elseif exp_mu eq 0;
        
        mu_ij	= zgama;
        
    endif;
    
    retp(mu_ij);
endp;

proc 1 = _IC(logL,N,k);
    local BIC, AIC, CAIC, HQIC;
    
    BIC   = (-2*logL) + (ln(N)*k);
    AIC   = (-2*logL) + (2*k);
    CAIC  = (-2*logL) + ((1+ln(N))*k);
    HQIC  = (-2*logL) + ((2*ln(ln(N)))*k);
    
    retp(BIC~AIC~CAIC~HQIC);
endp;


/* PROCEDURES FOR GRADIENTS */

proc 1 = _grad_beta_OP(mu_xbeta,d_ij,p_ij,x);
    local grad_beta, grad_betaJ;
    
    grad_beta		= pdfn(mu_xbeta);
    grad_betaJ		= pdfn(mu_xbeta[.,capJ-1] .* -1) .* -1;
    grad_beta		= grad_beta[.,1]~(grad_beta[.,2:cols(grad_beta)] - grad_beta[.,1:cols(grad_beta)-1])~grad_betaJ;
    grad_beta		= sumr(d_ij .* grad_beta);
    grad_beta		= (grad_beta .*. ones(1,cols(x))) .* (-x);
    grad_beta		= grad_beta ./ (sumr(d_ij .* p_ij) .*. ones(1,cols(x)));
    
    retp(grad_beta);
endp;



//proc 1 = _grad_beta_inf(d_ij,p_ij,p_ij_OP,X_inf,xb_inf,xb_inf_struct);
//    local grad_b_inf, p_aug0, p_aug1;
    
//    grad_b_inf	= (pdfn(xb_inf) .*. ones(1,capJ));
//    grad_b_inf 	= sumr(grad_b_inf .* p_ij_OP .* d_ij);
//    grad_b_inf	= (grad_b_inf .*. ones(1,cols(X_inf))) .* X_inf;
    
//    p_aug0  	= pdfn(-xb_inf) .* cdfn(-xb_inf_struct) .* d_ij[.,index0+1];
//    p_aug1  	= pdfn(-xb_inf) .* cdfn(xb_inf_struct) .* d_ij[.,index1+1];
    
//    p_aug0  	= (p_aug0 .*. ones(1,cols(X_inf))) .* -X_inf;
//    p_aug1  	= (p_aug1 .*. ones(1,cols(X_inf))) .* -X_inf;
    
//    grad_b_inf	= grad_b_inf + p_aug0;
//    grad_b_inf	= grad_b_inf + p_aug1;
    
//    p_ij		= sumr(p_ij .* d_ij) .*. ones(1,cols(X_inf));
    
//    grad_b_inf	= grad_b_inf ./ p_ij;
    
//    retp(grad_b_inf);
//endp;


//proc 1 = _grad_beta_inf_struct(d_ij,p_ij,X_struct,xb_inf,xb_inf_struct);
//    local grad_b_inf_struct, p_aug0, p_aug1;
    
//    p_aug0  	= cdfn(-xb_inf) .* pdfn(-xb_inf_struct) .* d_ij[.,index0+1] ./ p_ij[.,index0+1];
//    p_aug1  	= cdfn(-xb_inf) .* pdfn(xb_inf_struct) .* d_ij[.,index1+1] ./ p_ij[.,index1+1];
    
//    p_aug0  	= (p_aug0 .*. ones(1,cols(X_struct))) .* -X_struct;
//    p_aug1  	= (p_aug1 .*. ones(1,cols(X_struct))) .* X_struct;
    
//    grad_b_inf_struct	= p_aug0 + p_aug1;
    
//    retp(grad_b_inf_struct);
//endp;


//proc 1 = _grad_beta_TOP1(mu_xbeta,d_ij,p_ij,x,p_temp1,p_temp2);
//    local grad_beta1, grad_beta2, grad_beta3a, grad_beta3b, grad_beta3, grad_beta4a, grad_beta4b, grad_beta4, grad_betaJ, grad_beta;
    
//    grad_beta1		= (pdfn(mu_xbeta[.,1]) .* d_ij[.,1] ./ p_ij[.,1]);
//    grad_beta1		= (grad_beta1 .*. ones(1,cols(x))) .* (-x);
    
//    grad_beta2		= ((pdfn(mu_xbeta[.,2]) - pdfn(mu_xbeta[.,1])) .* d_ij[.,2] ./ p_ij[.,2] .* p_temp1[.,1]);
//    grad_beta2		= (grad_beta2 .*. ones(1,cols(x))) .* (-x);
    
//    grad_beta3a		= (pdfn(mu_xbeta[.,3]) - pdfn(mu_xbeta[.,2])) .* d_ij[.,3];
//    grad_beta3b		= ((pdfn(mu_xbeta[.,2]) - pdfn(mu_xbeta[.,1])) .* d_ij[.,3] .* p_temp1[.,2]);
//    grad_beta3		= (grad_beta3a + grad_beta3b) ./ p_ij[.,3];
//    grad_beta3		= (grad_beta3 .*. ones(1,cols(x))) .* (-x);
    
//    grad_beta4a		= (pdfn(mu_xbeta[.,4]) - pdfn(mu_xbeta[.,3])) .* d_ij[.,4];
//    grad_beta4a		= (grad_beta4a .*. ones(1,cols(x))) .* (-x);
//    grad_beta4b		= pdfn(-mu_xbeta[.,4]) .* d_ij[.,4] .* p_temp2[.,2];
//    grad_beta4b		= (grad_beta4b .*. ones(1,cols(x))) .* x;
//    grad_beta4		= (grad_beta4a + grad_beta4b) ./ (p_ij[.,4] .*. ones(1,cols(x)));
    
//    grad_betaJ		= pdfn(-mu_xbeta[.,capJ-1]) .* d_ij[.,5] ./ p_ij[.,5] .* p_temp2[.,1];
//    grad_betaJ		= (grad_betaJ .*. ones(1,cols(x))) .* x;
    
//    grad_beta		= grad_beta1 + grad_beta2 + grad_beta3 + grad_beta4 + grad_betaJ;
    
//    retp(grad_beta);
//endp;

//proc 1 = _grad_beta_TOP2(mu_xbeta,d_ij,p_ij,x,p_temp1,p_temp2);
//    local grad_beta1, grad_beta2, grad_beta3a, grad_beta3b, grad_beta3c, grad_beta3;
//    local grad_beta4a, grad_beta4b, grad_beta4c, grad_beta4, grad_betaJ, grad_beta;
    
//    grad_beta1		= (pdfn(mu_xbeta[.,1]) .* d_ij[.,1] ./ p_ij[.,1]);
//    grad_beta1		= (grad_beta1 .*. ones(1,cols(x))) .* (-x);
    
//    grad_beta2		= ((pdfn(mu_xbeta[.,2]) - pdfn(mu_xbeta[.,1])) .* d_ij[.,2] ./ p_ij[.,2] .* p_temp1[.,1]);
//    grad_beta2		= (grad_beta2 .*. ones(1,cols(x))) .* (-x);
    
//    grad_beta3a		= (pdfn(mu_xbeta[.,3]) - pdfn(mu_xbeta[.,2])) .* d_ij[.,3];
//    grad_beta3b		= ((pdfn(mu_xbeta[.,2]) - pdfn(mu_xbeta[.,1])) .* d_ij[.,3] .* p_temp1[.,2]);
//    grad_beta3c		= pdfn(mu_xbeta[.,4]) .* d_ij[.,3] .* p_temp2[.,1] .* -1;
//    grad_beta3		= (grad_beta3a + grad_beta3b + grad_beta3c) ./ p_ij[.,3];
//    grad_beta3		= (grad_beta3 .*. ones(1,cols(x))) .* (-x);
    
//    grad_beta4a		= (pdfn(mu_xbeta[.,4]) - pdfn(mu_xbeta[.,3])) .* d_ij[.,4];
//    grad_beta4a		= (grad_beta4a .*. ones(1,cols(x))) .* (-x);
//    grad_beta4b		= ((pdfn(mu_xbeta[.,2]) - pdfn(mu_xbeta[.,1])) .* d_ij[.,4] .* p_temp1[.,3]);
//    grad_beta4b		= (grad_beta4b .*. ones(1,cols(x))) .* (-x);
//    grad_beta4c		= pdfn(mu_xbeta[.,4]) .* d_ij[.,4] .* p_temp2[.,2];
//    grad_beta4c		= (grad_beta4c .*. ones(1,cols(x))) .* x;
//    grad_beta4		= (grad_beta4a + grad_beta4b + grad_beta4c) ./ (p_ij[.,4] .*. ones(1,cols(x)));
    
//    grad_betaJ		= pdfn(mu_xbeta[.,capJ-1]) .* d_ij[.,5] ./ p_ij[.,5] .* p_temp2[.,3];
//    grad_betaJ		= (grad_betaJ .*. ones(1,cols(x))) .* x;
    
//    grad_beta		= grad_beta1 + grad_beta2 + grad_beta3 + grad_beta4 + grad_betaJ;
    
//    retp(grad_beta);
//endp;


proc 1 = _grad_mu_OP(mu_xbeta,d_ij,p_ij);
    local grad_mu1, grad_mu2, grad_mu3, grad_mu4, grad_mu;
    
    grad_mu1		= (pdfn(mu_xbeta[.,1]) .* d_ij[.,1]) ./ p_ij[.,1];
    grad_mu1		= grad_mu1 + ((-pdfn(mu_xbeta[.,1]) .* d_ij[.,2]) ./ p_ij[.,2]);
    
    grad_mu2		= (pdfn(mu_xbeta[.,2]) .* d_ij[.,2]) ./ p_ij[.,2];
    grad_mu2		= grad_mu2 + ((-pdfn(mu_xbeta[.,2]) .* d_ij[.,3]) ./ p_ij[.,3]);
    
    if capJ eq 5;
        
        grad_mu3		= (pdfn(mu_xbeta[.,3]) .* d_ij[.,3]) ./ p_ij[.,3];
        grad_mu3		= grad_mu3 + ((-pdfn(mu_xbeta[.,3]) .* d_ij[.,4]) ./ p_ij[.,4]);
        
        grad_mu4		= (pdfn(mu_xbeta[.,4]) .* d_ij[.,4]) ./ p_ij[.,4];
        grad_mu4		= grad_mu4 + ((-pdfn(mu_xbeta[.,4] .* -1) .* d_ij[.,5]) ./ p_ij[.,5]);
        
        grad_mu			= grad_mu1~grad_mu2~grad_mu3~grad_mu4;
        
    elseif capJ eq 4;
        
        grad_mu3		= (pdfn(mu_xbeta[.,3]) .* d_ij[.,3]) ./ p_ij[.,3];
        grad_mu3		= grad_mu3 + ((-pdfn(mu_xbeta[.,3]) .* d_ij[.,4]) ./ p_ij[.,4]);
        
        grad_mu			= grad_mu1~grad_mu2~grad_mu3;
        
    elseif capJ eq 3;
        
        grad_mu			= grad_mu1~grad_mu2;
        
    endif;
    
    retp(grad_mu);
endp;


proc 1 = _grad_mu_TOP1(mu_xbeta,d_ij,p_ij,p_temp1,p_temp2);
    local grad_mu1, grad_mu2, grad_mu3, grad_mu4, grad_mu;
    
    grad_mu1		= (pdfn(mu_xbeta[.,1]) .* d_ij[.,1]) ./ p_ij[.,1];
    grad_mu1		= grad_mu1 + (-pdfn(mu_xbeta[.,1]) .* p_temp1[.,1] .* d_ij[.,2]) ./ p_ij[.,2];
    grad_mu1		= grad_mu1 + (-pdfn(mu_xbeta[.,1]) .* p_temp1[.,2] .* d_ij[.,3]) ./ p_ij[.,3];
    
    grad_mu2		= (pdfn(mu_xbeta[.,2]) .* p_temp1[.,1] .* d_ij[.,2]) ./ p_ij[.,2];
    grad_mu2		= grad_mu2 + (pdfn(mu_xbeta[.,2]) .* (p_temp1[.,2] - 1) .* d_ij[.,3]) ./ p_ij[.,3];
    
    grad_mu3		= (pdfn(mu_xbeta[.,3]) .* d_ij[.,3]) ./ p_ij[.,3];
    grad_mu3		= grad_mu3 + (-pdfn(mu_xbeta[.,3]) .* d_ij[.,4]) ./ p_ij[.,4];
    
    grad_mu4		= (pdfn(mu_xbeta[.,4]) .* (1 - p_temp2[.,2]) .* d_ij[.,4]) ./ p_ij[.,4];
    grad_mu4		= grad_mu4 + (-pdfn(mu_xbeta[.,4]) .* p_temp2[.,1] .* d_ij[.,5]) ./ p_ij[.,5];
    
    grad_mu			= grad_mu1~grad_mu2~grad_mu3~grad_mu4;
    
    retp(grad_mu);
endp;

proc 1 = _grad_mu_TOP2(mu_xbeta,d_ij,p_ij,p_temp1,p_temp2);
    local grad_mu1, grad_mu2, grad_mu3, grad_mu4, grad_mu;
    
    grad_mu1		= (pdfn(mu_xbeta[.,1]) .* d_ij[.,1]) ./ p_ij[.,1];
    grad_mu1		= grad_mu1 + (-pdfn(mu_xbeta[.,1]) .* p_temp1[.,1] .* d_ij[.,2]) ./ p_ij[.,2];
    grad_mu1		= grad_mu1 + (-pdfn(mu_xbeta[.,1]) .* p_temp1[.,2] .* d_ij[.,3]) ./ p_ij[.,3];
    grad_mu1		= grad_mu1 + (-pdfn(mu_xbeta[.,1]) .* p_temp1[.,3] .* d_ij[.,4]) ./ p_ij[.,4];
    
    grad_mu2		= (pdfn(mu_xbeta[.,2]) .* p_temp1[.,1] .* d_ij[.,2]) ./ p_ij[.,2];
    grad_mu2		= grad_mu2 + (((-pdfn(mu_xbeta[.,2])) + (pdfn(mu_xbeta[.,2]) .* p_temp1[.,2])) .* d_ij[.,3]) ./ p_ij[.,3];
    grad_mu2		= grad_mu2 + ((pdfn(mu_xbeta[.,2]) .* p_temp1[.,3] .* d_ij[.,4]) ./ p_ij[.,4]);
    
    grad_mu3		= (pdfn(mu_xbeta[.,3]) .* d_ij[.,3]) ./ p_ij[.,3];
    grad_mu3		= grad_mu3 + (-pdfn(mu_xbeta[.,3]) .* d_ij[.,4]) ./ p_ij[.,4];
    
    grad_mu4		= ((-pdfn(mu_xbeta[.,4])) .* p_temp2[.,1] .* d_ij[.,3]) ./ p_ij[.,3];
    grad_mu4		= grad_mu4 + ((pdfn(mu_xbeta[.,4]) + ((-pdfn(mu_xbeta[.,4])) .* p_temp2[.,2])) .* d_ij[.,4]) ./ p_ij[.,4];
    grad_mu4		= grad_mu4 + (((-pdfn(mu_xbeta[.,4])) .* p_temp2[.,3]) .* d_ij[.,5]) ./ p_ij[.,5];
    
    grad_mu			= grad_mu1~grad_mu2~grad_mu3~grad_mu4;
    
    retp(grad_mu);
endp;

proc 1 = _grad_b_tempered_TOP1(p_OP,d_ij,p_ij,xb_temp,X_inf);
    local grad_temp1, grad_temp2, grad_b_tempered;
    
    grad_temp1 = (pdfn(xb_temp) .* p_OP .* d_ij[.,1] ./ p_ij[.,1]) .*. ones(1,cols(X_inf));
    grad_temp1 = grad_temp1 .* X_inf;
    
    grad_temp2 = (pdfn(xb_temp) .* p_OP .* d_ij[.,2] ./ p_ij[.,2]) .*. ones(1,cols(X_inf));
    grad_temp2 = grad_temp2 .* (-X_inf);
    
    grad_b_tempered = grad_temp1 + grad_temp2;
    
    retp(grad_b_tempered);
endp;

proc 1 = _grad_b_tempered_TOP2(p_OP,d_ij,p_ij,xb_temp,X_inf);
    local grad_temp1, grad_temp2, grad_temp3, grad_b_tempered;
    
    grad_temp1 = (pdfn(xb_temp[.,1]) .* p_OP .* d_ij[.,1] ./ p_ij[.,1]) .*. ones(1,cols(X_inf));
    grad_temp1 = grad_temp1 .* (-X_inf);
    
    grad_temp2 = ((pdfn(xb_temp[.,2]) - pdfn(xb_temp[.,1])) .* p_OP .* d_ij[.,2] ./ p_ij[.,2]) .*. ones(1,cols(X_inf));
    grad_temp2 = grad_temp2 .* (-X_inf);
    
    grad_temp3 = (pdfn(xb_temp[.,2]) .* p_OP .* d_ij[.,3] ./ p_ij[.,3]) .*. ones(1,cols(X_inf));
    grad_temp3 = grad_temp3 .* X_inf;
    
    grad_b_tempered = grad_temp1 + grad_temp2 + grad_temp3;
    
    retp(grad_b_tempered);
endp;


proc 1 = _grad_tempered_mu1_TOP2(p_OP,d_ij,p_ij,xb_temp);
    local grad_temp1, grad_temp2, grad_temp3, grad_mu1_tempered;
    
    grad_temp1 = (pdfn(xb_temp[.,1]) .* p_OP .* d_ij[.,1] ./ p_ij[.,1]);
    grad_temp2 = ((pdfn(xb_temp[.,2]) - pdfn(xb_temp[.,1])) .* p_OP .* d_ij[.,2] ./ p_ij[.,2]);
    grad_temp3 = (-pdfn(xb_temp[.,2]) .* p_OP .* d_ij[.,3] ./ p_ij[.,3]);
    
    grad_mu1_tempered = grad_temp1 + grad_temp2 + grad_temp3;
    
    retp(grad_mu1_tempered);
endp;

proc 1 = _grad_tempered_mu2_TOP2(p_OP,d_ij,p_ij,xb_temp,exp_mu_2);
    local grad_temp1, grad_temp2, grad_mu2_tempered;
    
    grad_temp1 = (pdfn(xb_temp) .* exp_mu_2 .* p_OP .* d_ij[.,1] ./ p_ij[.,1]);
    grad_temp2 = (pdfn(-xb_temp) .* (-exp_mu_2) .* p_OP .* d_ij[.,2] ./ p_ij[.,2]);
    
    grad_mu2_tempered = grad_temp1 + grad_temp2;
    
    retp(grad_mu2_tempered);
endp;

proc 1 = _grad_V_HOPIT(sig_V_inv,mu_xbetaV,p_ijV,V_ij);
    local grad_V_con_1, grad_V_temp, grad_V_con_J, grad_V_con;
    
    grad_V_con_1	= -sig_V_inv .* pdfn(mu_xbetaV[.,1]);
    grad_V_temp		= -sig_V_inv .* (pdfn(mu_xbetaV[.,2:capJ-1]) - pdfn(mu_xbetaV[.,1:capJ-2]));
    grad_V_con_J	= sig_V_inv .* pdfn(mu_xbetaV[.,capJ-1] .* -1);
    grad_V_con		= grad_V_con_1~grad_V_temp~grad_V_con_J;
    grad_V_con		= sumr((grad_V_con ./ p_ijV) .* V_ij);
    
    retp(grad_V_con);
endp;

proc 1 = _grad_theta_HOPIT(mu_xbetaV,sig_V_inv,p_ijV,V_ij);
    local grad_theta_V_1, grad_theta_temp, grad_theta_V_J, grad_theta_V;
    
    grad_theta_V_1	= (mu_xbetaV[.,1] ./ sig_V_inv) .* pdfn(mu_xbetaV[.,1]);
    grad_theta_temp	= ((mu_xbetaV[.,2:capJ-1] ./ sig_V_inv) .* pdfn(mu_xbetaV[.,2:capJ-1])) -
        ((mu_xbetaV[.,1:capJ-2] ./ sig_V_inv) .* pdfn(mu_xbetaV[.,1:capJ-2]));
    grad_theta_V_J	= ((mu_xbetaV[.,capJ-1] ./ sig_V_inv) .* -1) .* pdfn(mu_xbetaV[.,capJ-1] .* -1);
    grad_theta_V	= grad_theta_V_1~grad_theta_temp~grad_theta_V_J;
    grad_theta_V	= sumr((grad_theta_V ./ p_ijV) .* V_ij);
    
    retp(grad_theta_V);
endp;


proc 1 = _grad_g_j(gama,z,mu_xbeta,p_ij,d_ij,sig_V_inv);
    local j1, j_mid, J, j_all, grad_g, grad_gV1, grad_gV2, grad_gV3, sel, j_rep;
    
    j1			= pdfn(mu_xbeta[.,1]);
    j_mid		= pdfn(mu_xbeta[.,2:cols(mu_xbeta)]) - pdfn(mu_xbeta[.,1:cols(mu_xbeta)-1]);
    J			= (pdfn(mu_xbeta[.,cols(mu_xbeta)] .* -1)) .* -1;
    j_all		= (j1~j_mid~J) .* sig_V_inv;
    
    if exp_mu eq 1;
        
        j_all		= (j_all .* (exp(z*gama) .*. ones(1,cols(p_ij)))) ./ p_ij;
        grad_g		= sumr(j_all .* d_ij);
        grad_g		= (grad_g .*. ones(1,cols(z))) .* z;
        
    elseif exp_mu eq 0;
        
        j_all		= (j_all ./ p_ij) .* d_ij;
        grad_g		= sumr(j_all) .* z;
        
    endif;
    
    retp(grad_g);
endp;

proc 1 = _grad_g_j_linear(Z,mu_xbeta,p_ij,d_ij,sig_V_inv);
    local j, j_mid, j_all, grad_g, grad_gV1, grad_gV2, grad_gV3, sel, j_rep;
    
    j_all		= pdfn(mu_xbeta) .* sig_V_inv;
    j_all		= sumr((j_all ./ p_ij) .* d_ij);
    grad_g		= sumr(j_all) .* Z;
    
    retp(grad_g);
endp;

proc 1 = _grad_g_j_linearJ(x,mu_xbeta,p_ij,d_ij,sig_V_inv);
    local j, j_mid, j_all, grad_g, grad_gV1, grad_gV2, grad_gV3, sel, j_rep;
    
    //j_all		= ((pdfn(mu_xbeta .* -1)) .* -1) .* sig_V_inv;
    j_all		= pdfn(mu_xbeta .* -1) .* sig_V_inv;
    //j_all		= pdfn(mu_xbeta) .* sig_V_inv;
    j_all		= sumr((j_all ./ p_ij) .* d_ij);
    grad_g		= sumr(j_all) .* x;
    
    retp(grad_g);
endp;

proc 1 = _grad_g1_TOP1(gama,Z,mu_xbeta,p_ij,d_ij,p_temp1,p_temp2);
    local j1, j2, j3a, j3b, j3, j4a, j4b, j4, j_mid, J, j_all, grad_g;
    
    j1			= pdfn(mu_xbeta[.,1])  .* exp(Z*gama) .* d_ij[.,1] ./ p_ij[.,1];
    j1			= (j1 .*. ones(1,cols(Z))) .* Z;
    
    j_mid		= (pdfn(mu_xbeta[.,2:cols(mu_xbeta)]) - pdfn(mu_xbeta[.,1:cols(mu_xbeta)-1]));
    
    j2			= j_mid[.,1] .* p_temp1[.,1] .* exp(Z*gama) .* d_ij[.,2] ./ p_ij[.,2];
    j2			= (j2 .*. ones(1,cols(Z))) .* Z;
    
    j3a			= ((j_mid[.,1] .* p_temp1[.,2] .* exp(Z*gama)) .*. ones(1,cols(Z))) .* Z;
    j3b			= ((j_mid[.,2]  .* exp(Z*gama)) .*. ones(1,cols(Z))) .* Z;
    j3			= (j3a + j3b) ./ (p_ij[.,3] .*. ones(1,cols(Z)));
    j3			= j3 .* (d_ij[.,3] .*. ones(1,cols(Z)));
    
    j4a			= ((j_mid[.,3]  .* exp(Z*gama)) .*. ones(1,cols(Z))) .* Z;
    j4b			= ((pdfn(mu_xbeta[.,4]) .* p_temp2[.,2]  .* exp(Z*gama)) .*. ones(1,cols(Z))) .* (-Z);
    j4			= (j4a + j4b) ./ (p_ij[.,4] .*. ones(1,cols(Z)));
    j4			= j4 .* (d_ij[.,4] .*. ones(1,cols(Z)));
    
    J			= (pdfn(mu_xbeta[.,4]) .* p_temp2[.,1]  .* exp(Z*gama) ./ p_ij[.,5]) .* (-Z);
    J			= J .* (d_ij[.,5] .*. ones(1,cols(Z)));
    
    grad_g		= j1 + j2 + j3 + j4 + J;
    
    retp(grad_g);
endp;


/*
proc 1 = _grad_g1_TOP1(Z,mu_xbeta,p_ij,d_ij,p_temp1,p_temp2);
local j1, j2, j3a, j3b, j3, j4a, j4b, j4, j_mid, J, j_all, grad_g;

j1			= pdfn(mu_xbeta[.,1]) .* d_ij[.,1] ./ p_ij[.,1];
j1			= (j1 .*. ones(1,cols(Z))) .* Z;

j_mid		= pdfn(mu_xbeta[.,2:cols(mu_xbeta)]) - pdfn(mu_xbeta[.,1:cols(mu_xbeta)-1]);

j2			= j_mid[.,1] .* p_temp1[.,1] .* d_ij[.,2] ./ p_ij[.,2];
j2			= (j2 .*. ones(1,cols(Z))) .* Z;

j3a			= ((j_mid[.,1] .* p_temp1[.,2]) .*. ones(1,cols(Z))) .* Z;
j3b			= (j_mid[.,2] .*. ones(1,cols(Z))) .* Z;
j3			= (j3a + j3b) ./ (p_ij[.,3] .*. ones(1,cols(Z)));
j3			= j3 .* (d_ij[.,3] .*. ones(1,cols(Z)));

j4a			= (j_mid[.,3] .*. ones(1,cols(Z))) .* Z;
j4b			= ((pdfn(mu_xbeta[.,4]) .* p_temp2[.,2]) .*. ones(1,cols(Z))) .* (-Z);
j4			= (j4a + j4b) ./ (p_ij[.,4] .*. ones(1,cols(Z)));
j4			= j4 .* (d_ij[.,4] .*. ones(1,cols(Z)));

J			= (pdfn(mu_xbeta[.,4]) .* p_temp2[.,1] ./ p_ij[.,5]) .* (-Z);
J			= J .* (d_ij[.,5] .*. ones(1,cols(Z)));

grad_g		= j1 + j2 + j3 + j4 + J;

retp(grad_g);
endp;
*/

/*
proc 1 = _grad_g1_TOP2(Z,mu_xbeta,p_ij,d_ij,p_temp1,p_temp2);
local j1, j2, j3a, j3b, j3c, j3, j4a, j4b, j4c, j4, j_mid, J, j_all, grad_g;

j1			= (pdfn(mu_xbeta[.,1]) ./ p_ij[.,1]) .* d_ij[.,1];
j1			= (j1 .*. ones(1,cols(Z))) .* Z;

j2			= ((pdfn(mu_xbeta[.,2]) - pdfn(mu_xbeta[.,1])) .* p_temp1[.,1] ./ p_ij[.,2]) .* d_ij[.,2];
j2			= (j2 .*. ones(1,cols(Z))) .* Z;

j3a			= ((pdfn(mu_xbeta[.,3]) - pdfn(mu_xbeta[.,2])) .*. ones(1,cols(Z))) .* Z;
j3b			= (((pdfn(mu_xbeta[.,2]) - pdfn(mu_xbeta[.,1])) .* p_temp1[.,2]) .*. ones(1,cols(Z))) .* Z;
j3c			= ((pdfn(-mu_xbeta[.,cols(mu_xbeta)]) .* p_temp2[.,1]) .*. ones(1,cols(Z))) .* (-Z);
j3			= (j3a + j3b + j3c) ./ (p_ij[.,3] .*. ones(1,cols(Z)));
j3			= j3 .* (d_ij[.,3] .*. ones(1,cols(Z)));

j4a			= ((pdfn(mu_xbeta[.,4]) - pdfn(mu_xbeta[.,3])) .*. ones(1,cols(Z))) .* Z;
j4b			= (((pdfn(mu_xbeta[.,2]) - pdfn(mu_xbeta[.,1])) .* p_temp1[.,3]) .*. ones(1,cols(Z))) .* Z;
j4c			= ((pdfn(mu_xbeta[.,4]) .* p_temp2[.,2]) .*. ones(1,cols(Z))) .* (-Z);
j4			= (j4a + j4b + j4c) ./ (p_ij[.,4] .*. ones(1,cols(Z)));
j4			= j4 .* (d_ij[.,4] .*. ones(1,cols(Z)));

J			= (pdfn(-mu_xbeta[.,4]) .* p_temp2[.,3] ./ p_ij[.,5]) .* (-Z);
J			= J .* (d_ij[.,5] .*. ones(1,cols(Z)));

grad_g		= j1 + j2 + j3 + j4 + J;

retp(grad_g);
endp;
*/

proc 1 = _grad_g1_TOP2(gama,Z,mu_xbeta,p_ij,d_ij,p_temp1,p_temp2);
    local j1, j2, j3a, j3b, j3c, j3, j4a, j4b, j4c, j4, j_mid, J, j_all, grad_g;
    
    j1			= (pdfn(mu_xbeta[.,1]) .* exp(Z*gama) ./ p_ij[.,1]) .* d_ij[.,1];
    j1			= (j1 .*. ones(1,cols(Z))) .* Z;
    
    j2			= ((pdfn(mu_xbeta[.,2]) - pdfn(mu_xbeta[.,1])) .* exp(Z*gama) .* p_temp1[.,1] ./ p_ij[.,2]) .* d_ij[.,2];
    j2			= (j2 .*. ones(1,cols(Z))) .* Z;
    
    j3a			= (((pdfn(mu_xbeta[.,3]) - pdfn(mu_xbeta[.,2])) .* exp(Z*gama)) .*. ones(1,cols(Z))) .* Z;
    j3b			= (((pdfn(mu_xbeta[.,2]) - pdfn(mu_xbeta[.,1])) .* exp(Z*gama) .* p_temp1[.,2]) .*. ones(1,cols(Z))) .* Z;
    j3c			= ((pdfn(-mu_xbeta[.,cols(mu_xbeta)]) .* exp(Z*gama) .* p_temp2[.,1]) .*. ones(1,cols(Z))) .* (-Z);
    j3			= (j3a + j3b + j3c) ./ (p_ij[.,3] .*. ones(1,cols(Z)));
    j3			= j3 .* (d_ij[.,3] .*. ones(1,cols(Z)));
    
    j4a			= (((pdfn(mu_xbeta[.,4]) - pdfn(mu_xbeta[.,3])) .* exp(Z*gama)) .*. ones(1,cols(Z))) .* Z;
    j4b			= (((pdfn(mu_xbeta[.,2]) - pdfn(mu_xbeta[.,1])) .* exp(Z*gama) .* p_temp1[.,3]) .*. ones(1,cols(Z))) .* Z;
    j4c			= (((pdfn(mu_xbeta[.,4]) .* p_temp2[.,2]) .* exp(Z*gama)) .*. ones(1,cols(Z))) .* (-Z);
    j4			= (j4a + j4b + j4c) ./ (p_ij[.,4] .*. ones(1,cols(Z)));
    j4			= j4 .* (d_ij[.,4] .*. ones(1,cols(Z)));
    
    J			= (pdfn(-mu_xbeta[.,4]) .* exp(Z*gama) .* p_temp2[.,3] ./ p_ij[.,5]) .* (-Z);
    J			= J .* (d_ij[.,5] .*. ones(1,cols(Z)));
    
    grad_g		= j1 + j2 + j3 + j4 + J;
    
    retp(grad_g);
endp;

proc 1 = _grad_g2_TOP1(gama,Z,mu_xbeta,p_ij,d_ij,p_temp1,p_temp2);
    local j2, j3a, j3b, j3, j4a, j4b, j4, j_mid, J, j_all, grad_g;
    
    j2			= pdfn(mu_xbeta[.,2]) .* exp(Z*gama) .* p_temp1[.,1] .* d_ij[.,2] ./ p_ij[.,2];
    j2			= (j2 .*. ones(1,cols(Z))) .* Z;
    
    j_mid		= pdfn(mu_xbeta[.,2:cols(mu_xbeta)]) - pdfn(mu_xbeta[.,1:cols(mu_xbeta)-1]);
    
    j3a			= ((pdfn(mu_xbeta[.,2]) .* exp(Z*gama) .* p_temp1[.,2]) .*. ones(1,cols(Z))) .* Z;
    j3b			= ((j_mid[.,2] .* exp(Z*gama)) .*. ones(1,cols(Z))) .* Z;
    j3			= (j3a + j3b) ./ (p_ij[.,3] .*. ones(1,cols(Z)));
    j3			= j3 .* (d_ij[.,3] .*. ones(1,cols(Z)));
    
    j4a			= ((pdfn(mu_xbeta[.,4]) .* exp(Z*gama) .* p_temp2[.,2]) .*. ones(1,cols(Z))) .* (-Z);
    j4b			= ((j_mid[.,3] .* exp(Z*gama)) .*. ones(1,cols(Z))) .* Z;
    j4			= (j4a + j4b) ./ (p_ij[.,4] .*. ones(1,cols(Z)));
    j4			= j4 .* (d_ij[.,4] .*. ones(1,cols(Z)));
    
    J			= (pdfn(mu_xbeta[.,4]) .* exp(Z*gama) .* p_temp2[.,1] ./ p_ij[.,5]) .* (-Z);
    J			= J .* (d_ij[.,5] .*. ones(1,cols(Z)));
    
    grad_g		= j2 + j3 + j4 + J;
    
    retp(grad_g);
endp;

proc 1 = _grad_g2_TOP2(gama,Z,mu_xbeta,p_ij,d_ij,p_temp1,p_temp2);
    local j2, j3a, j3b, j3c, j3, j4a, j4b, j4c, j4, j_mid, J, j_all, grad_g;
    
    j2			= pdfn(mu_xbeta[.,2]) .* exp(Z*gama) .* p_temp1[.,1] .* d_ij[.,2] ./ p_ij[.,2];
    j2			= (j2 .*. ones(1,cols(Z))) .* Z;
    
    j3a			= (((pdfn(mu_xbeta[.,3]) - pdfn(mu_xbeta[.,2])) .* exp(Z*gama)) .*. ones(1,cols(Z))) .* Z;
    j3b			= ((pdfn(mu_xbeta[.,2]) .* exp(Z*gama) .* p_temp1[.,2]) .*. ones(1,cols(Z))) .* Z;
    j3c			= ((pdfn(mu_xbeta[.,4]) .* exp(Z*gama) .* p_temp2[.,1]) .*. ones(1,cols(Z))) .* (-Z);
    j3			= (j3a + j3b +j3c) ./ (p_ij[.,3] .*. ones(1,cols(Z)));
    j3			= j3 .* (d_ij[.,3] .*. ones(1,cols(Z)));
    
    j4a			= (((pdfn(mu_xbeta[.,4]) - pdfn(mu_xbeta[.,3])) .* exp(Z*gama)) .*. ones(1,cols(Z))) .* Z;
    j4b			= ((pdfn(mu_xbeta[.,2]) .* exp(Z*gama) .* p_temp1[.,3]) .*. ones(1,cols(Z))) .* Z;
    j4c			= ((pdfn(mu_xbeta[.,4]) .* exp(Z*gama) .* p_temp2[.,2]) .*. ones(1,cols(Z))) .* (-Z);
    j4			= (j4a + j4b + j4c) ./ (p_ij[.,4] .*. ones(1,cols(Z)));
    j4			= j4 .* (d_ij[.,4] .*. ones(1,cols(Z)));
    
    J			= ((pdfn(mu_xbeta[.,4]) .* exp(Z*gama) .* p_temp2[.,3]) ./ p_ij[.,5]) .* (-Z);
    J			= J .* (d_ij[.,5] .*. ones(1,cols(Z)));
    
    grad_g		= j2 + j3 + j4 + J;
    
    
    retp(grad_g);
endp;

proc 1 = _grad_g3_TOP1(gama,Z,mu_xbeta,p_ij,d_ij,p_temp2);
    local j3a, j3b, j3, j4a, j4b, j4, j_mid, J, j_all, grad_g;
    
    j3			= pdfn(mu_xbeta[.,3]) .* exp(Z*gama) .* d_ij[.,3] ./ p_ij[.,3];
    j3			= (j3 .*. ones(1,cols(Z))) .* Z;
    
    j_mid		= pdfn(mu_xbeta[.,2:cols(mu_xbeta)]) - pdfn(mu_xbeta[.,1:cols(mu_xbeta)-1]);
    
    j4a			= ((pdfn(mu_xbeta[.,4]) .* exp(Z*gama) .* p_temp2[.,2]) .*. ones(1,cols(Z))) .* (-Z);
    j4b			= ((j_mid[.,3] .* exp(Z*gama)) .*. ones(1,cols(Z))) .* Z;
    j4			= (j4a + j4b) ./ (p_ij[.,4] .*. ones(1,cols(Z)));
    j4			= j4 .* (d_ij[.,4] .*. ones(1,cols(Z)));
    
    J			= (pdfn(mu_xbeta[.,4]) .* exp(Z*gama) .* p_temp2[.,1] ./ p_ij[.,5]) .* (-Z);
    J			= J .* (d_ij[.,5] .*. ones(1,cols(Z)));
    
    grad_g		= j3 + j4 + J;
    
    retp(grad_g);
endp;

proc 1 = _grad_g3_TOP2(gama,Z,mu_xbeta,p_ij,d_ij,p_temp2);
    local j3a, j3b, j3, j4a, j4b, j4, j_mid, J, j_all, grad_g;
    
    j3a			= pdfn(mu_xbeta[., 3]) .* exp(Z*gama);
    j3a			= (j3a .*. ones(1,cols(Z))) .* Z;
    j3b			= ((pdfn(mu_xbeta[., 4]) .* exp(Z*gama) .* p_temp2[., 1]) .*. ones(1, cols(Z))) .* (-Z);
    j3			= (j3a + j3b) ./ (p_ij[., 3] .*. ones(1, cols(Z)));
    j3			= j3 .* (d_ij[., 3] .*. ones(1, cols(Z)));
    
    j_mid		= pdfn(mu_xbeta[., 2:cols(mu_xbeta)]) - pdfn(mu_xbeta[., 1:cols(mu_xbeta)-1]);
    
    j4a			= ((pdfn(mu_xbeta[., 4]) .* exp(Z*gama) .* p_temp2[., 2]) .*. ones(1, cols(Z))) .* (-Z);
    j4b			= ((j_mid[., 3] .* exp(Z*gama)) .*. ones(1, cols(Z))) .* Z;
    j4			= (j4a + j4b) ./ (p_ij[., 4] .*. ones(1, cols(Z)));
    j4			= j4 .* (d_ij[., 4] .*. ones(1, cols(Z)));
    
    J			= (pdfn(mu_xbeta[., 4]) .* exp(Z*gama) .* p_temp2[., 3] ./ p_ij[., 5]) .* (-Z);
    J			= J .* (d_ij[., 5] .*. ones(1, cols(Z)));
    
    grad_g		= j3 + j4 + J;
    
    retp(grad_g);
endp;

proc 1 = _grad_g4_TOP1(gama, Z, mu_xbeta, p_ij, d_ij, p_temp2);
    local j4a, j4b, j4, J, j_all, grad_g;
    
    j4a			= ((pdfn(mu_xbeta[., 4]) .* exp(Z*gama) .* p_temp2[., 2]) .*. ones(1,cols(Z))) .* (-Z);
    j4b			= ((pdfn(mu_xbeta[., 4]) .* exp(Z*gama)) .*. ones(1,cols(Z))) .* Z;
    j4			= (j4a + j4b) ./ (p_ij[., 4] .*. ones(1,cols(Z)));
    j4			= j4 .* (d_ij[., 4] .*. ones(1,cols(Z)));
    
    J			= (pdfn(mu_xbeta[., 4]) .* exp(Z*gama) .* p_temp2[., 1] ./ p_ij[., 5]) .* (-Z);
    J			= J .* (d_ij[., 5] .*. ones(1,cols(Z)));
    
    grad_g		= j4 + J;
    
    retp(grad_g);
endp;

proc 1 = _grad_g4_TOP2(gama, Z, mu_xbeta, p_ij, d_ij, p_temp2);
    local j3, j4a, j4b, j4, J, j_all, grad_g;
    
    j3			= ((pdfn(mu_xbeta[., 4]) .* exp(Z*gama) .* p_temp2[., 1]) .*. ones(1, cols(Z))) .* (-Z);
    j3			= j3 ./ (p_ij[., 3] .*. ones(1, cols(Z)));
    j3			= j3 .* (d_ij[., 3] .*. ones(1, cols(Z)));
    
    j4a			= ((pdfn(mu_xbeta[., 4]) .* exp(Z*gama) .* p_temp2[., 2]) .*. ones(1, cols(Z))) .* (-Z);
    j4b			= ((pdfn(mu_xbeta[., 4]) .* exp(Z*gama)) .*. ones(1, cols(Z))) .* Z;
    j4			= (j4a + j4b) ./ (p_ij[., 4] .*. ones(1, cols(Z)));
    j4			= j4 .* (d_ij[., 4] .*. ones(1, cols(Z)));
    
    J			= (pdfn(mu_xbeta[., 4]) .* exp(Z*gama) .* p_temp2[., 3] ./ p_ij[., 5]) .* (-Z);
    J			= J .* (d_ij[., 5] .*. ones(1, cols(Z)));
    
    grad_g		= j3 + j4 + J;
    
    retp(grad_g);
endp;

proc 1 = _grad_gJ(gama,x,mu_xbeta,p_ij,d_ij,sig_V_inv);
    local j1, j_all, J, grad_g;
    
    j1			= pdfn(mu_xbeta);
    J			= (pdfn(mu_xbeta .* -1)) .* -1;
    j_all		= (j1~J) .* sig_V_inv;
    j_all		= (j_all .* (exp(x*gama) .*. ones(1,cols(p_ij)))) ./ p_ij;
    
    grad_g		= sumr(j_all .* d_ij);
    grad_g		= (grad_g .*. ones(1,cols(x))) .* x;
    
    retp(grad_g);
endp;

