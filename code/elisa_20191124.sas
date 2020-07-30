
%let indir =U:\Consulting\KEL\Fox\McEntee\input;
%let outdir =U:\Consulting\KEL\Fox\McEntee\output;
libname elisa "&indir." access=readonly;
libname out "&outdir.";
filename mylst "&outdir.\my_listing.lst";
filename mylog "&outdir.\my_listing.log";



data elisa0; set elisa.elisa ; 
dead = (dead_alive = 0);
run;

data elisa1; set elisa.elisa_20191124; run;

data elisa; merge elisa0 elisa1; by cat; run;



%let cont_dv_iris = len_sag pelvis_transverse Pelvis_saggital proxuret_cm parench_cm_trans parench_cm_sag
Ratio_PS_RS_sag Ratio_PS_RS_trans;
%let cat_iv_iris = pre_IRIS dis_IRIS short_IRIS long_IRIS;

%let cont_vars = age weight duration_CS pre_BUN pre_creat pre_SDMA 
pre_K pre_PCV pre_USG Dis_BUN Dis_creat Dis_SDMA Delta_op_BUN 
Delta_op_Creat Delta_op_SDMA Short_BUN Short_Creat Short_SDMA 
Long_BUN Long_creat Long_SDMA Hist_BUN Hist_creat Hist_SDMA 
Delta_BUN Delta_creat Delta_SDMA FU_pelvis_cm Follow_up_time_months;
%let cat_vars = sex left_right Long_SUB_flush_patency dead 
Cause_death Death_renal Death_mayberenal Death_nonrenal SUB_related;

/* hist_SDMA - No data */
%let cont_iv = Age Weight Left_Right pre_USG pre_PCV len_sag pelvis_transverse Pelvis_saggital proxuret_cm parench_cm_trans
parench_cm_sag Ratio_PS_RS_sag Ratio_PS_RS_trans Duration_CS hist_BUN hist_creat  pre_BUN pre_Creat pre_SDMA;
%let cont_dv = pre_BUN pre_creat pre_SDMA pre_PCV short_BUN short_creat short_SDMA delta_BUN
delta_creat long_BUN long_creat long_SDMA delta_op_BUN delta_op_creat delta_op_SDMA parench_cm_trans parench_cm_sag
pre_IRIS dis_IRIS short_IRIS long_IRIS;
%let cat_iv = Left_Right margin_irreg renal_asym stricture_stone cause_ID ips_nephro cysto contra_nephro;
%let cat_dv = dead death_renal margin_irreg renal_asym;
;


*** Continuous *************;


%macro means_desc(ds);
%let vars = &cont_vars.;
%let word_cnt=%sysfunc(countw(&vars));
%do i = 1 %to &word_cnt.;
proc means data =  &ds. noprint;
var %scan(&vars.,&i);

output out= mean mean(%scan(&vars.,&i))=;
output out= median median(%scan(&vars.,&i))=;
output out= std std(%scan(&vars.,&i))= ;
output out= min min(%scan(&vars.,&i))= ;
output out= max max(%scan(&vars.,&i))=;
run;
data mean (rename = (%scan(&vars.,&i) = mean _freq_ = n)); set mean (drop =  _type_);
data median (rename = %scan(&vars.,&i) = median); set median (keep = %scan(&vars.,&i)); 
data std (rename = %scan(&vars.,&i) = std); set std (keep = %scan(&vars.,&i)); 
data min (rename = %scan(&vars.,&i) = min); set min (keep = %scan(&vars.,&i)); 
data max (rename = %scan(&vars.,&i) = max); set max (keep = %scan(&vars.,&i)); 
data all_mean_%scan(&vars.,&i); length ds $30.; merge mean median std min max; 
ds = "%scan(&vars.,&i)"; run;
data all_mean; set all_mean all_mean_%scan(&vars.,&i); run;
%end;
	%mend;

data all_mean; set _null_; run;
%means_desc(elisa);



*** Categorical *************;

* Descriptive Categorical One-Way Freq ; 
%macro freqs1(ds);
%let vars = &cat_vars;
%let word_cnt=%sysfunc(countw(&vars));
%do i = 1 %to &word_cnt.;
proc freq data =  &ds.;
tables %scan(&vars.,&i) / chisq;
ods output onewayfreqs = owf OneWayChiSq= owc;
run;
data freq (rename = f_%scan(&vars.,&i) = var_mid); 
 set owf (drop = %scan(&vars.,&i)); run;
data chi (drop = name1 label1); set owc; where name1 = "P_PCHI"; run;
data freq_%scan(&vars.,&i); length var_mid cValue1 $50.; merge freq chi; by table; 
data freq_tot; set freq_tot freq_%scan(&vars.,&i); var = left(var_mid); run;
%end;
	%mend;

	data freq_tot; set _null_; run;
	%freqs1(elisa);


/* Regression */

%macro ols_mac(ds);
%let count_iv = %sysfunc(countw(&cont_iv.));
%let count_dv = %sysfunc(countw(&cont_dv.));
%do i = 1 %to &count_iv.; * &count_iv.;
%do j = 1 %to &count_dv.; * &count_dv. ;
%put IV:&count_iv. DV:&count_dv.;
%put IV_s:"%scan(&cont_dv,&j.)" DV_s:"%scan(&cont_iv,&j.)";
%if "%scan(&cont_dv,&j.)" ne "%scan(&cont_iv,&i.)" %then %do;

/*%end;*/
/*	%end;*/
/*%mend;*/
/*	%ols_mac;*/
proc reg data=&ds.;
  model %scan(&cont_dv,&j.) = %scan(&cont_iv,&i.) ;
  output out=res rstudent=r /*h=lev cookd=cd dffits=dffit*/ ;
  ods output ParameterEstimates=pe_out FitStatistics = fs_out (keep=label2 nvalue2);
run;
quit;
data fs (drop=label2 nvalue2); set fs_out; where label2 = "R-Square"; rho = sqrt(nvalue2); r2 = nvalue2; run;
data pe_1;
length variable iv dv $50.;
set pe_out;
iv = "%scan(&cont_iv,&i.)";
dv = "%scan(&cont_dv,&j.)";
where variable not in ("Intercept");
sig = (probt lt 0.05); 
run;


data _null_; set fs; call symput("rho",rho); run;
proc power; onecorr sides = 1 alpha = 0.05 corr = &rho. power = 0.8 ntotal = .; ods output output = power_out; run;
%put &rho.;

data power_out (drop= variable Analysis Index Sides NominalAlpha NominalPower NullCorr NPartialVars Error Info);
length variable iv dv $50.; set power_out; iv = "%scan(&cont_iv,&i.)"; dv = "%scan(&cont_dv,&j.)";  run;


data fs_pe; merge pe_1 power_out;
data pe;
set pe fs_pe;
run;


data _null_;
   set pe_out;
   if _n_ = 1 then call symput('Int', put(estimate, BEST6.));    
   else            call symput('Slope', put(estimate, BEST6.));  
run;

proc sgplot data=&ds.;
   title "Regression Line with Slope and Intercept";
   reg y=%scan(&cont_dv,&j.) x=%scan(&cont_iv,&i.);
   inset "Intercept = &Int" "Slope = &Slope" / 
         border title="Parameter Estimates" position=topleft;
run;

* Examine residuals for normality ;
proc univariate data=res plots plotsize=30 normal;
  var %scan(&cont_iv,&i.);
  ods output TestsForNormality=norm_1;
run;
data norm_1;
length varname iv dv $50.;
set norm_1;
iv = "%scan(&cont_iv,&i.)";
dv = "%scan(&cont_dv,&j.)";
run;
data norm;
set norm norm_1;
run;


%end;
	%end;
		%end;
%mend;

proc printto print=mylst log=mylog new; run;


data norm; set _null_; data pe; set _null_; data pe_norm; set _null_; data npar; set _null_; run;
%ols_mac(elisa);




/* ANOVA */

%macro anova(ds);
%let count_iv = %sysfunc(countw(&cat_iv.));
%let count_dv = %sysfunc(countw(&cont_dv.));
%do i = 1 %to &count_iv.; * &count_iv.;
%do j = 1 %to &count_dv.; * &count_dv.;
%put IV:&count_iv. DV:&count_dv.;
%put IV_s:"%scan(&cont_dv,&j.)" DV_s:"%scan(&cont_iv,&j.)";
%put &count_iv. &count_dv.;
proc mixed data=&ds.;
   class %scan(&cat_iv,&i.) ;
   model %scan(&cont_dv,&j.) = %scan(&cat_iv,&i.) / solution;
   lsmeans %scan(&cat_iv,&i.) / pdiff=all ; 
   ods output lsmeans=lsmeans_x (drop = df tvalue probt) 
	diffs = diffs_x (keep = estimate probt %scan(&cat_iv,&i.)  _%scan(&cat_iv,&i.)); run;
%if %sysfunc(exist(diffs_x)) %then %do;
data diffs (rename = (estimate = Compr_mean)); 
set diffs_x;  mergeme = 1; 
run;
data lsmeans (rename = (estimate = lsmean effect = iv)); 
length dv /*iv*/ $50.;
set lsmeans_x  ; 
/*iv = "%scan(&cat_iv,&i.)";*/
dv = "%scan(&cont_dv,&j.)";
lsmean_level = left(trim(put(%scan(&cat_iv,&i.),5.0)));
mergeme = 1;
run;
data prep (drop=mergeme); length iv $50.; merge diffs (drop =  %scan(&cat_iv,&i.) _%scan(&cat_iv,&i.)) 
lsmeans (drop =  %scan(&cat_iv,&i.)); by mergeme; run;
data allaov ; retain dv iv lsmean_level lsmean stderr; 
set allaov prep; run;
	%end;
/*	proc datasets lib=work nolist nodetails; delete diffs_x lsmeans_x diffs lsmeans ; run;*/
/* Get normality */
proc reg data=&ds.;
  model %scan(&cont_dv,&j.) = %scan(&cat_iv,&i.) ;
  output out=res rstudent=r /*h=lev cookd=cd dffits=dffit*/ ;
  ods output ParameterEstimates=pe_1;
run;
quit;
data pe_1;
length variable iv dv $50.;
set pe_1;
iv = "%scan(&cat_iv,&i.)";
dv = "%scan(&cont_dv,&j.)";
where variable not in ("Intercept");
run;
data pe;
set pe pe_1;
run;
* Examine residuals for normality ;
proc univariate data=res plots plotsize=30 normal;
  var %scan(&cat_iv,&i.);
  ods output TestsForNormality=norm_1;
run;
data norm_1;
length varname iv dv $50.;
set norm_1;
iv = "%scan(&cat_iv,&i.)";
dv = "%scan(&cont_dv,&j.)";
run;
data norm;
set norm norm_1;
run;

%end;
	%end;
%mend;

data allaov; set _null_; data lsmeans; set _null_; data diffs; set _null_; run;
data norm; set _null_; data pe; set _null_; data pe_norm; set _null_; data npar; set _null_; run;

%anova(elisa);

/* IRIS ANOVA */


%macro anova_iris(ds,iv);
/*%let cat_iv = &cat_iv_iris.;*/
%let cont_dv = &cont_dv_iris.;
/*%let count_iv = %sysfunc(countw(&cat_iv.));*/
%let count_dv = %sysfunc(countw(&cont_dv.));
/*%do i = 1 %to &count_iv.; * &count_iv.;*/
%do j = 1 %to &count_dv.; * &count_dv.;
/*%put IV:&count_iv. DV:&count_dv.;*/
/*%put IV_s:"%scan(&cont_dv,&j.)" DV_s:"%scan(&cont_iv,&j.)";*/
/*%put &count_iv. &count_dv.;*/
proc mixed data=&ds.;
/*   class %scan(&cat_iv,&i.) ;*/
class &iv.;
/*   model %scan(&cont_dv,&j.) = %scan(&cat_iv,&i.) / solution;*/
/*lsmeans %scan(&cat_iv,&i.) / pdiff=all ; */
/*model %scan(&cont_dv,&j.) = %scan(&cat_iv,&i.) / solution;*/
model %scan(&cont_dv,&j.) = &iv. ;
lsmeans &iv. / pdiff=all ; 
   ods output lsmeans=lsmeans_x (drop = df tvalue probt) 
/*	diffs = diffs_x (keep = estimate probt %scan(&cat_iv,&i.)  _%scan(&cat_iv,&i.)); run;*/
	diffs = diffs_x (keep = estimate probt &iv. _&iv.); run;
%if %sysfunc(exist(diffs_x)) %then %do;
data diffs (rename = (estimate = Compr_mean)); 
set diffs_x;  mergeme = 1; 
/*lsmean_level1 = left(trim(put(%scan(&cat_iv,&i.),5.0)));*/
/*lsmean_level2 = left(trim(put(%scan(_&cat_iv,&i.),5.0)));*/
lsmean_level1 = left(trim(put((&iv.),5.0)));
lsmean_level2 = left(trim(put((_&iv.),5.0)));
run;
data lsmeans (rename = (estimate = lsmean effect = iv)); 
length dv /*iv*/ $50.;
set lsmeans_x  ; 
/*iv = "%scan(&cat_iv,&i.)";*/
dv = "%scan(&cont_dv,&j.)";
/*lsmean_level = left(trim(put(%scan(&cat_iv,&i.),5.0)));*/
lsmean_level = left(trim(put((&iv.),5.0)));
mergeme = 1;
run;
/*data prep (drop=mergeme); length iv $50.; merge diffs (drop =  %scan(&cat_iv,&i.) _%scan(&cat_iv,&i.)) */
/*lsmeans (drop =  %scan(&cat_iv,&i.)); by mergeme; run;*/
data prep (drop=mergeme); length iv $50.; merge diffs (drop =  &iv. _&iv.) 
lsmeans (drop =  &iv.); by mergeme; run;
data allaov ; retain dv iv lsmean_level lsmean stderr; 
set allaov prep; run;
	%end;
	%end;
%mend;

data allaov; set _null_; data lsmeans; set _null_; data diffs; set _null_; run;
/*data norm; set _null_; data pe; set _null_; data pe_norm; set _null_; data npar; set _null_; run;*/

%anova_iris(elisa,pre_IRIS);
%anova_iris(elisa,dis_IRIS);
%anova_iris(elisa,short_IRIS);
%anova_iris(elisa,long_IRIS);





/* Logistic Continuous */



%macro log_cont(ds,dv);
%let count_iv = %sysfunc(countw(&cont_iv.));
%do i = 1 %to &count_iv.; * &count_iv.;
/*ods trace on;*/
proc logistic data = &ds.;
/*class Left_Right (ref = '1') ;*/
model &dv. (event = '1') = %scan(&cont_iv,&i.)  / RISKLIMITS /* LACKFIT RSQUARE PARMLABEL*/;
ods output ParameterEstimates=pe Association=a (keep=label2 nvalue2) CLoddsWald=or (drop=effect);
/*ods trace off;*/
%put "&dv." "%scan(&cont_iv,&i.)";
data a2(keep=c_stat); set a (rename=nvalue2=c_stat); if label2 = "c"; run;
data pe2 (rename=(estimate=log_estimate ProbChiSq=p_value)); length variable $50.; set pe (drop=df); if variable = "%scan(&cont_iv,&i.)"; run;
data or_a_pe; length dv $50.; merge or a2 pe2; dv = "&dv.";
data or_a_pe; retain dv Variable OddsRatioEst LowerCL UpperCL p_value c_stat Unit log_est stdErr WaldChiSq _ESTTYPE_; set or_a_pe;  run;
data or_all; set or_all or_a_pe; run;
run;
	%end;
%mend;

data or_all; set _null_; run;
%log_cont(elisa,dead);
%log_cont(elisa,death_renal);
%log_cont(elisa,margin_irreg);
%log_cont(elisa,renal_asym);


%macro log_cat(ds,dv);
%let count_iv = %sysfunc(countw(&cat_iv.));
%do i = 1 %to &count_iv.; * &count_iv.;
%if "&dv." ne "%scan(&cat_iv,&i.)" %then %do;
proc logistic data = &ds.;
class %scan(&cat_iv,&i.) (ref = '0') ;
model &dv. (event = '1') = %scan(&cat_iv,&i.)  / RISKLIMITS /* LACKFIT RSQUARE PARMLABEL*/;
ods output ParameterEstimates=pe Association=a (keep=label2 nvalue2) CLoddsWald=or (drop=effect);
/*ods trace off;*/
%put "&dv." "%scan(&cat_iv,&i.)";
data a2(keep=c_stat); set a (rename=nvalue2=c_stat); if label2 = "c"; run;
data pe2 (rename=(estimate=log_estimate ProbChiSq=p_value)); length variable $50.; set pe (drop=df); if variable = "%scan(&cat_iv,&i.)"; run;
data or_a_pe; length dv $50.; merge or a2 pe2; dv = "&dv.";
data or_a_pe; retain dv Variable OddsRatioEst LowerCL UpperCL p_value c_stat Unit log_est stdErr WaldChiSq _ESTTYPE_; set or_a_pe;  run;
data or_all; set or_all or_a_pe; run;
run;
	%end;
		%end;
%mend;

data or_all; set _null_; run;
%log_cat(elisa,dead);
%log_cat(elisa,death_renal);
%log_cat(elisa,margin_irreg);
%log_cat(elisa,renal_asym);


/* Two Way Freqs */


%macro chi(ds,analysis,r_q_n,research,iv1);
%let catvars =  &cat_iv.;
%let count=%sysfunc(countw(&catvars.));
%do i = 1 %to &count.;
%if "&iv1." ne "%scan(&cat_iv,&i.)" %then %do;
proc freq data =  &ds.;
tables &iv1.*%scan(&catvars,&i.) / chisq;
ods output CrossTabFreqs = freqs ChiSq= chi FishersExact = fish ;run;
run;
data freqs_in; set freqs (drop =  _table_ rowpercent colpercent); 
where _type_ in ("11" "00"); drop _type_; run;


data ctf_&research. (drop =  &iv1. %scan(&catvars,&i.)); retain r_q_n analysis research_q; length r_q_n research_q $50.; set freqs_in; 
analysis = "&analysis."; research_q = "&research."; r_q_n = "&r_q_n."; 
iv1 = put(&iv1.,5.0); iv2 = put(%scan(&catvars,&i.),5.0);
run;
data ctf_all; set ctf_all ctf_&research.; run;

	%if %sysfunc(exist(fish)) %then %do;
data f_&research. (rename=nvalue1=P_FISH); set fish (keep =  table name1 nvalue1); 
where name1 in ("P_TABLE"); drop name1; run;
data f_p; set f_p f_&research.; run;
	%end;
	%if %sysfunc(exist(chi)) %then %do;
data chi_&research. (rename=prob=P_CHI); set chi (keep =  table statistic prob); 
where statistic in ("Chi-Square");  drop statistic; 
r_q_n = "&r_q_n."; run;
data chi_p; set chi_p chi_&research.; run;
	%end;

proc datasets library=work nolist nodetails; delete chi fish freq p_&research. chi_&research.; 
run;
%end;
	%end;
	%mend;
	data ctf_all; set _null_; data f_p; set _null_; data chi_p; set _null_; run;
%chi(elisa,chi,999,placeholder,dead);
%chi(elisa,chi,999,placeholder,death_renal);
%chi(elisa,chi,999,placeholder,margin_irreg);
%chi(elisa,chi,999,placeholder,renal_asym);

