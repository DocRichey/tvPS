# tvPS
Repository for time-varying propensity score code

/* Welcome, here I share SAS code and the data structure used in the comparison of time-varying propensity scores (tvPS) with exact
matching approaches - manuscript at: https://bmcmedresmethodol.biomedcentral.com/counter/pdf/10.1186/s12874-024-02391-3.pdf 

Any questions please feel free to contact me at MorganRichey@gmail.com - I'll do my best to respond quickly*/

/* Introduction:  If this is your first time working with PHREG (or just need a refresher), please review this excellent guide to Cox PH regression from UCLA stats: https://stats.oarc.ucla.edu/sas/seminars/sas-survival/ .  

If possible, I HIGHLY recommend downloading the example dataset and running through this entire tutorial if time permits */

  /* Data structure:  This approach utilizes the counting process approach, which requires a data structure that includes a start/stop covariate for each time point.  While this approach may appear to require more work than the programming statement approach, I found that QC was simpler with this method.  See section 2.1, paragraph 2 (and second table) in https://stats.oarc.ucla.edu/sas/seminars/sas-survival/     /* 

/* Step 0:  Decide on which covariates should be included in your propensity score model, I recommend the below paper by Brookhart et al that concludes that confounders and predictors of the outcome should be included (they do not recommend including covariates that solely predict the exposure)
https://pmc.ncbi.nlm.nih.gov/articles/PMC1513192/pdf/nihms8838.pdf
*/

/* Step 1:  Estimating a time-varying propensity score in a cohort of cases and eligible controls, ties handled with Efron method */

/* Example from tvPS paper - categorical covariates go in class statement, continuous covariates in the model statement */
proc phreg data = ready outest=haz ;
class 
diabetic (ref="0")
case (ref="0")
SEX (ref="0")
Insulin (ref="0")
DEPR (ref="0")
DYSLIPID (ref="0")
EATINGDO (ref="0")
FATTYLIVER (ref="0")
MI (ref="0")
/param=ref;
model (start, stop)*case(0) = diabetic SEX Insulin DEPR DYSLIPID FATTYLIVER GERD MI Rbmi Rgagne RVTERISK Rage 
/ TIES = EFRON RL; 
 BASELINE OUT=SURVS SURVIVAL=S CUMHAZ=CumHaz;;
 id pid; /* note this is the individual ID covariate - all individuals should have a unique ID (pid=person ID) for all encounters - even if they change status (e.g. act as a control for half the study, then later act as a case) */
 run;

/* Carefully review the parameter estimates for each covariate, see UCLA PHREG review and associated literature for diagnostics, be wary of large hazard estimates - I highly recommend comparing crosstabs between cases and controls for covariates that have large parameter estimates.  Verify that covariates that are excellent predictors of case or control status are both internally and externally valid. */

/* Step 2:  Assign hazard values to cases and covariates by using output from PHREG -0  */


data puts;
set parms;
if _n_ =1 then call symput('diabetic', estimate);
if _n_ =2 then call symput('sex', estimate);
if _n_ =3 then call symput('insulin', estimate);
if _n_ =4 then call symput('DEPR', estimate);
if _n_ =5 then call symput('DYSLIPID', estimate);
if _n_ =6 then call symput('FATTYLIVER', estimate);
if _n_ =7 then call symput('MI', estimate);
if _n_ =8 then call symput('Rbmi', estimate);
if _n_ =9 then call symput('Rgagne', estimate);
if _n_ =10 then call symput('rVTERISK', estimate);
if _n_ =11 then call symput('Rage', estimate);
run;


 /* Assign hazard values to each participant in the dataset from which the propensity score was estimated*/
 
data hazards5;
set ready;
hazdiabetic= diabetic*&diabetic;
hazsex= sex * &sex;
hazinsulin= insulin * &insulin;
hazDEPR= depr* &DEPR;
hazDYSLIPID= dyslipid * &DYSLIPID;
hazFATTYLIVER= fattyliver* &FATTYLIVER;
hazMI= MI * &MI;
hazRbmi= rbmi * &Rbmi;
hazRgagne= rgagne * &Rgagne;
hazrVTERISK= rvterisk * &rVTERISK;
hazRage= rage * &Rage;
haz=(hazdiabetic+hazsex+hazinsulin+hazDEPR+hazDYSLIPID+hazFATTYLIVER+hazMI+hazRbmi+hazRgagne+hazrVTERISK+hazRage);
if haz<0 then haz=.00001; /* We must constrain the hazard score to between 0 and 1, so the lowest hazard score becomes a propensity score of .00001 */
run;

proc sort data=hazardsa5;
by case;
run;

/* Identify the largest hazard score - this will become a propensity score of 1.0 */
proc means data=hazards5 ;
output out=hazmax max=hazmax;
var haz;
run;

/* Define the maximum hazard */
data _null;
set hazmax;
call symput('hazmax',hazmax);
run;

data hazscore5;
set hazards5;
hazscore=(haz/&hazmax);
run;

/* Match cases to controls 1:1, with replacement - 
1.  note that your case identifier must go in the class statement
2.  You can force exact matching on some covariates, in our case we required a control to have a healthcare encounter to be within a 6-month window prior to a case's index surgical encounter
3.  If you are doing 1:1 matching, I have observed SAS matching a single case to > 1 control if the control propensity scores are (a) close enough and (b) controls are sufficiently plentiful.  As a QC check, I recommend double check the number of participants with a _MatchID and ensure it is double the number of cases.    
* /

proc psmatch data=hazscore5;
class case matchtime;
psdata treatvar=case(Treated="1") ps=hazscore;
match method=replace(k=1) distance=PS stat=lps caliper=0.25 exact=(matchtime);
assess ps var=(hazscore BMI Male Diabetes insulin Depression Dyslipidemia FattyLiver MI rgagne rVTErisk stroke VTE GERD HiatalHernia AFIB CIRRHOSIS EATINGDisorder);
output out(obs=match)=haz100replace lps=_Lps matchid=_MatchID;
run;

/* 
