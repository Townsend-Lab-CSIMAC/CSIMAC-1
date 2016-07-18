#if !defined(PRFCluster_H)
#define  PRFCluster_H


#define CRI(a) ((a)==0?"BIC":(a)==1?"AIC":"AICc") //criteria, models

//Declare several points of description of the program
#define NAME "CSI-MAC"	//Program name
#define FULLNAME "Cancer Selection Intensity - Model Averaged Clustering"
#define VERSION    "1.0"
#define LASTUPDATE "July 18th, 2016"
#define FUNCTION "Estimate selection intensity for each single site in coding sequences by using divergence data in tumors."
#define REFERENCE "Zi-Ming Zhao, Sachith Gullapalli, Ning Li and Jeffrey P. Townsend. (2016)"



#define IR_H 50      //Initial value of r
#define IR_L -50    // Initial r
#define ER 0.0001   //f(r)<=ER, this r is the one we are trying to find. Exact r.
#define MinDx 0.4   //f(r)<=MinDx, use the optimal gamma meeting this criteria even though ER is not reached
#define JMAX 1000  //Maximum number of iterations of Newton-Raphson's method to approximate gamma
#define GAMMALIST 3801 //Total number of gamma included in the lookup table for gamma from -100 to 100, with 0.1 interval and 0.01 from -10 to 10; totally 3803 per line, with one n=12 and one pmuc.


#define N_Random 300000 //Choose models randomly, how many models are chosen. No effects to the speed and memory to the program by having 1000 instead of 10000.


#include "base.h"
#include <map>
#include <time.h>
#ifdef __linux__
#include <sys/sysinfo.h>
#endif

using namespace std;

//Four different models are used for calculation, namely InL, AIC, AICc, and BIC. All models are different ways of clustering for a particular site
class CandidateModels {
 public:
  float CW; //criterion and weight under the specific criterion
  double AICweight; //AIC weight
  double LogLikelihood;
  

  long pos_start, pos_end; // start and end position of the sequence
  float p0, pc; // probability of the non-cluster and cluster
  long cs, ce; // cluster start position and end position
  
  double InL0, InL; // LogLikelihood
  double AIC0, AIC; // AIC
  double AICc0, AICc; // AICc
  double BIC0, BIC; // BIC
  


  /***************************************************
  * Function: For all possible candidate models for clustering
  * Input Parameter:criterion weight - CW, criteria type, AIC, BIC, or AICc; sequence start and end positions; cluster start and end positions; probability of non-cluster and cluster regions
  * Output:
  * Return Value:
  ***************************************************/

  CandidateModels (double AICweight, float CW, long pos_start, long pos_end, long cs, long ce, float p0, float pc, double LogLikelihood){
	  this->AICweight = AICweight;
	  this->CW = CW;
    this->pos_start = pos_start;		this->pos_end = pos_end;
    this->cs = cs;		this->ce = ce;
    this->p0 = p0;		this->pc = pc;
    this->LogLikelihood = LogLikelihood;
  }


  /***************************************************
  * Function: For best-fit selected models for clustering
  * Input Parameter:
  * Output:
  * Return Value:
  ***************************************************/
  CandidateModels (long pos_start, long pos_end, long cs, long ce, float p0, float pc, double InL0, double InL,double AIC0, double AIC, double AICc0, double AICc, double BIC0, double BIC) {
    this->pos_start = pos_start;		this->pos_end = pos_end;
    this->cs = cs;		this->ce = ce;
    this->p0 = p0;		this->pc = pc;
    this->InL0 = InL0;		this->InL = InL;
    this->AIC0 = AIC0;		this->AIC = AIC;
    this->AICc0 = AICc0;		this->AICc = AICc;
    this->BIC0 = BIC0;		this->BIC = BIC;
  }  
};
struct less_than_Models
{
    bool operator() (const CandidateModels& model1, const CandidateModels& model2)
    {
       if (model1.AICweight < model2.AICweight)
    	   return true;
       else return false;
    }
};


/***************************************************
* Function: class of lookup table contents with the recurrent number, lowerCIlambda, upperCIlambda
* Input Parameter:
* Output:
* Return Value:
***************************************************/
class CIRecurrentLookup {
 public:
  int count;
  double lowerCIlambda;
  double upperCIlambda;

  CIRecurrentLookup(int count, double lowerCIlambda, double upperCIlambda){
    this->count = count;
	this->lowerCIlambda = lowerCIlambda;
    this->upperCIlambda = upperCIlambda;
  }
};


/***************************************************
* Function: class of recurrent sites with site position and the recurrent number
* Input Parameter:
* Output:
* Return Value:
***************************************************/
class RecurrentSitesCounts {
 public:
  int sites;
  int counts;

  RecurrentSitesCounts(int sites, int counts){
    this->sites = sites;
    this->counts = counts;
  }
};

/***************************************************
* Function: structure containing weight and probability for each position
* Input Parameter:
* Output:
* Return Value:
***************************************************/
class CI {
 public:
  float weight;
  float p;
  float LogLikelihood;

  CI(float weight, float p, float LogLikelihood){
    this->weight = weight;		
    this->p = p; // this p is the mutation rate of the site being in the cluster or not in the cluster
    this->LogLikelihood=LogLikelihood;
  }
};

/***************************************************
* Function: define a structure named SiteModels with a variable called sms; for a particular site, weight and probability are called using CI structure.
* Input Parameter:
* Output:
* Return Value:
***************************************************/

struct SiteModels{
  long pos;
  vector<CI> sms; //the name of the vector struct is sms, SiteModelSelection, it contains weight and probability.
};


struct more_than_CI
{
    bool operator() (const CI& struct1, const CI& struct2)
    {
       if ((struct1.weight > struct2.weight) || ((struct1.weight == struct2.weight) && (struct1.p > struct2.p)))
    	   return true;
       else return false;
    }
};


/***************************************************
* Function: structure defining gamma models including weight and gamma value for each site.
* Input Parameter:
* Output:
* Return Value:
***************************************************/
class rModels {
 public:
  float weight;
  float r;
  
  rModels(float weight, float r){
    this->weight=weight;
    this->r=r;
  }
};

/***************************************************
* Function:
* Input Parameter:
* Output:
* Return Value:
***************************************************/
class cPRFCluster: public Base {	
 public:	
  cPRFCluster();
  ~cPRFCluster();	
  
  //Main function
  int Run(int argc, const char*argv[]);

  //Lookup table from Mathematica to get the integration for gamma calculation
  //vector<double> GammaLookupTable(int tmp_n, string input_f_name, vector<double> Fn);
  int LambdaCILookupTable(string input_f_name);
  int GetRecurrentList(string input_f);
  //Clustering using maximum likelihood
  int RunML(vector<string> div_seq);
  

  //AIC/BIC based on sub sequence
  int ClusterSubSeq(int pos_start, int pos_end, char symbol='S', struct SiteModels *p=NULL);
  
  
  int init(long N);
  int ModelAveraging(long pos_start, long pos_end, long cs, long ce, double p0, double pc, double min_cri,struct SiteModels *p);
  int CI_MA(struct SiteModels *p,long N);

  int SitePRF(int species_n, long N);
  int rh_SitePRF(int species_n, double uhr, long N);
  int rc_SitePRF(int species_m, double ucr, long N);
  int recurrent_SiteGamma(int tumor_num, double ucr, int RecurrentNum, int Site, ostringstream* myout);
  int recurrent_SiteGammaCI(int tumor_num, double ucr, double CI, int Site, int upOrlow, ostringstream* myout);//upOrlow=0, low; upOrlow=1, up
  int rh_gx_SitePRF(int species_n, double uhr, long N);

  int SiteNI(long N);

  double DivergentTime(long species_n, long N);
  double PolymorphismSynonymousRate(long species_n, long N);
  double CancerSynonymousRate(long species_m, long N);
  double ReplacementRate(double ratio_NS, double syn_rate);


  //double ClusterPRF(long species_n, double pr, double dr);
	

  int CIr_stochastic(struct SiteModels *dr, long N);
  void CIr_stochastic_threaded(struct SiteModels *dr, long N, long i, time_t time_start1, std::ostringstream* myout);

  int BubbleSort(struct SiteModels *p,long site);
  long RandomModel_Num(struct SiteModels *p,long site);

  vector<double> rWeightSums;
  //vector<double> pWeightSums;
  vector<double> RandomModel_NumFastInit(struct SiteModels *p,long site);
  long RandomModel_NumFast(const vector<double> &pWeightSums);

  
  int CIr_exact(struct SiteModels *dr, long N);
  void CIr_exact_threaded(struct SiteModels *dr, long N, long i, std::ostringstream* myout);

  double CIs_PRF(double p_pr, double p_dr,long species_n);
  double CIs_rc_PRF(double ucr, double p_dr,long species_m);


  int CI_UpLow_rc(long site,double min_weight,vector<rModels> vec_rModels_c_indiv, ostringstream* myout);

  int output(long N);

 protected:
  
  double BinomialProb(long n, long i);
  double factorial(int n); 
  double LogLikelihoodCluster(long cs, long ce, long start, long end, char symbol='S');
  double LogLikelihoodNonCluster(long cs, long ce, long start, long end, char symbol='S');


  
  //int nw (from pos_start to pos_end, all of the nonsynonymous)
  //int nc (from cs to ce, all of the nonsynonymous)
  double getp0pc_MK(int pos_start, int pos_end, int cs, int ce, float &p0, float &pc, int nw, int nc);


  //flag_N_pol=0 means there is all A/T/G/C in the polymophism data, otherwise there is Ns or other ambiguous nucleotide
  //flag_N_div=0 means there is all A/T/G/C in the divergence data, otherwise there is Ns or other ambiguous nucleotide
  int flag_N_pol;
  int flag_N_div;
  
  double SilentRate;
  int SilentRate_flag;

  //Get Polymorphism Synonymous and Replacements
  string getPolSysRep(vector<string> seq);

  //Get Divergence Synonymous and Replacements
  string getDivSysRep(vector<string> pol, vector<string> div);
  
  //Consider the sites other than A, T, G or C 
  int ReplaceCodon4PolSys(vector <string>& codon, vector <string>& codon_other);
  int ReplaceCodon4DivSys(vector <string>& codon, vector <string>& codon_other, string& div_codon);
  char ReplaceSite(char codon_other_symbol, int pointer_site[], string symbol);

  //
  long getDifference(string seq, int pos_start, int pos_end, char symbol='1');
 
  
  //Parse parameters
  bool parseParameter(int argc, const char* argv[]);
  //Show help information for cMAC-PRF
  void showHelpInfo();
  
  
 public:
  //A vector of selected models
  vector<CandidateModels> vec_SelectedModels;
  
  //A vector of all candidate models
  vector<CandidateModels> vec_AllModels;
  
  //Store all the models from pr*dr 
  //Instead of p, it is estimated r.
  //Just for one site each time
  vector<rModels> vec_rModels;
  vector<rModels> vec_rModels_h; // for human polymorphism
  vector<rModels> vec_rModels_c; // for cancer
  
  vector <CIRecurrentLookup> LambdaCIs;//the class to hold the k, lowerCI_Lambda, UpperCI_Lambda


  //Lower and Upper CI for gamma
  vector<double> vec_lower_r;
  vector<double> vec_upper_r;
   
  //Lower and Upper CI for gamma in human polymorphism
  vector<double> vec_lower_r_h;
  vector<double> vec_upper_r_h;

  //Lower and Upper CI for gamma in cancer
  vector<double> vec_lower_r_c;
  vector<double> vec_upper_r_c;
  
  vector<double> vec_lower_r_c_r;
  vector<double> vec_upper_r_c_r;
  
  //rate by model selection
  vector<double> vec_MS_rate;
  //H=Hot Spot, C=Cold Spot, default=-
  vector<char> vec_spot;
  //rate by model averaging
  vector<double> vec_MA_rate;

  //Lower and Upper CI for MA rate
  vector<double> vec_lower_rate;
  vector<double> vec_upper_rate;
  
  //Vectors. Store the divergence time and selection coefficient r for each site.
  vector<double> vec_time;
  vector<double> vec_r;
  vector<double> vec_r_h; //vector of gamma for human polymorphism
  vector<double> vec_r_c; //vector of gamma for cancer
  vector<double> vec_r_c_r; //vector of gamma for cancer, recurrent only
  
  int TotalRecurrentCount;
  int TotalRecurrentSite;
  int TotalReplacementSite;
  int TotalSilentSite;
  double uhs;// synonymous human polymorphism germline mutation rate
  double ucs;// synonymous cancer divergence mutation rate

  double uhr; // replacement human polymorphism mutation rate
  double ucr; // replacement cancer divergence mutation rate

  
  //ps/pr/ds/dr>1 && flag_found==0, it means it could not reject the null model
  //It means when cs=0, ce=N-1, then cri==cri0, and cri could never be less than cri0 (0-N-1,para=0)
  int flag_found_pr;
  int flag_found_dr;
  int flag_found_ps;
  int flag_found_ds;
  
  vector<CandidateModels> vec_SelectedModels_ps;
  vector<CandidateModels> vec_SelectedModels_pr;
  vector<CandidateModels> vec_SelectedModels_ds;
  vector<CandidateModels> vec_SelectedModels_dr;
  
  vector<CandidateModels> vec_AllModels_ps;
  vector<CandidateModels> vec_AllModels_pr;
  vector<CandidateModels> vec_AllModels_ds;
  vector<CandidateModels> vec_AllModels_dr;
  
  // MS synonymous and replacements for polymorphism and divergence
  vector<double> vec_MS_rate_ps;
  vector<double> vec_MS_rate_pr;
  vector<double> vec_MS_rate_ds;
  vector<double> vec_MS_rate_dr;

  // MA synonymous and replacements for polymorphism and divergence
  vector<double> vec_MA_rate_ps;
  vector<double> vec_MA_rate_pr;
  vector<double> vec_MA_rate_ds;
  vector<double> vec_MA_rate_dr;

  // lower synonymous and replacements for polymorphism and divergence
  vector<double> vec_lower_rate_ps;
  vector<double> vec_lower_rate_pr;
  vector<double> vec_lower_rate_ds;
  vector<double> vec_lower_rate_dr;

  // upper synonymous and replacements for polymorphism and divergence
  vector<double> vec_upper_rate_ps;
  vector<double> vec_upper_rate_pr;
  vector<double> vec_upper_rate_ds;
  vector<double> vec_upper_rate_dr;

  //For estimating Neutrality Index
  vector<double> vec_NI;

 protected:

  //required parameters tumor sequences number and polymorphism sequences number in the format of string and int
  string tumor_num_s;
  string polymorphism_num_s;

  int tumor_num;
  int polymorphism_num;
  
  //polymorphism consensus file with Replacement and Sysnonymous labeled.
  string pol_cons_seqfile;
  vector<string> pol_cons_seqname;
  vector<string> pol_cons_seq;
  string pol_codon_consensus;

  //Cancer divergence consensus file with Replacement and Sysnonymous labeled.
  string div_cons_seqfile;
  vector<string> div_cons_seqname;
  vector<string> div_cons_seq;
  string div_codon_consensus;

  string recurrent_file;
  vector <RecurrentSitesCounts> Recurrents;

  int output_format_num;

  //Choose the genetic code for this species
  int genetic_code;

  //Criterion, 0=BIC, 1=AIC, 2=AICc, 3=LRT
  int criterion_type;

  //Confidence interval for model averaging, default=0.95
  float confidence_interval;
  
  //Threshold quantile for upper and lower bound of confidence interval for model averaging, default=0.025
  float quantile_for_CI;

  //Only model selection, default=0
  int MS_only;


  //Calculate 95% confidence intervals for Model Averaging ci_ma, default 0
  int ci_ma;
  
  //Estimate selection coefficient r, default 1
  int Do_r_estimate;
  
  //Calculate 95% confidence intervals for selection coefficient Do_ci_r, default 1
  int Do_ci_r;
  
  //Algorithm for calculating 95% confidence interval for selection coefficient Do_ci_r_exact, default 0
  int Do_ci_r_exact;
  
  //ZMZ 04/28/2016  added option - regional gamma only
  //regional only gamma or weighted gamma of regional and recurrent; default is weighted gamma
  int regional_gamma_only;
  
  //Cluster synonymous sites for polymorphism and divergence. 
  //Default is 0 without showing the clustering results of synonymous from polymorphism and divergent sequences
  //User could choose 1 to show the clustering results of synonymous.
  int Do_Synonymous_Cluster;
  
  //Scale of the gene, it should be 1,3,6,9,12,15...3*n, the true sequence length is the times of the scale and the given length.
  int Scale;

  //Divergent time; User could input their own species divergence time
  double Div_time;
  
  //Default is 1 which means to replace ambiguous nucleotide (N, R, Y, etc.) with the most frequently used nucleotide in other sequences
  //Otherwise is 0, see this codon as a gap
  int Nuc_replace;

  //Estimate the Neutrality Index for each site [integer, optional], {0: NOT estimate Neutrality Index || 1: estimate Neutrality Index}, default=0. 
  //Neutrality Index=(PR/DR)/(PS/DS)
  int NI_estimate;

 private:
};

#endif


