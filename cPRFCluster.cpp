/*
To change: multiple hits considered poisson rate

ZMZ 07/20/2016
Fixed bug: output Scale; 
if (output_format_num==1 || (output_format_num==0 && Scale>1)) //(output_format_num==0 and Scale==3) and (output_format_num==1 and Scale==1) are exactly the same output position
    {
		if (output_format_num==0 and Scale>1) {Scale=Scale/3;}//for amino acid, scale it down 3 times to make it comparable to nucleotide Scale
		for(long i=0; i<N; i++) {//each position after scaled
			for (long j=0;j<Scale;j++){ // JT: fix this Scale/3 = 1/3 - ZMZ fixed adjusting amino acid scales down 3 times to be comparable to nucleotide and expanding each position i to i*Scale+j (0<j<Scale).
				cout.width(width);cout<<Scale*i+j+1<<"\t"; //convert to the real amino acid or nucleotide position of the gene, if Scale==1, i; if Scale=6, i=10, new position would be 60-65. 
Fixed Seven compiling warnings:
changed from unassigned to int for int position1 = str.find("	");
|| and && parenthesis
deleted (k,tmp_l,tmp_u)=(0,0.0,0.0); use k=0; tmp_l=tmp_u=0.0;
Added output error messages and updated codes in parseParameters 

Changed p=0, r=0 in CIr_stochastic_threaded and CIs_rc_PRF: CIr_stochastic_threaded p=0, r=0; CIs_rc_PRF, dr=0, r=0.
Excluded one site cluster: 
In ClusterSubSeq, added ce-cs>1 in 'if (cri <= cri0 && ce-cs>1) {//add the condition to exclude the cluster of one site'
Modified: Used only subregional models to calculate gamma.
Modified: p=0, gamma=0; and remove conditions of vec_r_c[i*3]==0 or vec_r_c[i]==0, that gives NULL gamma

Added GeneLength public parameter, when n=0 and N=GeneLength, quit cluster; else keep all cluster models for sub-regions.
Modified: In LogLikelihoodNonCluster and LogLikelihoodCluster, return 0 when n=0 and lambda=0.
EachSiteModels: all models from the regions with clusters are kept, including models from divide and conquer ClusterSubSeq
//pointer[i].sms.clear();// silent pointer[i].sms.clear() in EachSiteModels, since cumulative models are used in ClusterSubSeq for each site.

ModelAveraging and MS_only==0:
Removed the option MS_only==0 as a condition in ClusterSubSeq, since site models are calculated anyway. 
Changed ModelAveraging to EachSiteModels to be more accurate and intuitive. 

//ZMZ 07/18/2016. 
//Modified: In LogLikelihoodNonCluster and LogLikelihoodCluster, to prevent log(0) case, get the pseudo-lambda: when n=0, lambda = (double)1/(N_cluster_ScaledBack+1);
Fixed bug: 
added symbol in LogLikelihoodNonCluster and LogLikelihoodCluster - bug in Synonymous count in likelihood calculation
LogLikelihoodNonCluster(long cs, long ce, long start, long end, char symbol)
LogLikelihoodCluster(long cs, long ce, long start, long end, char symbol)
Fixed bug in Recurrent count in LogLikelihoodNonCluster, from && to ||, (Recurrents[i].sites < cs || Recurrents[i].sites > ce)
Fixed bug in LogLikelihoodCluster and LogLikelihoodNonCluster: added double, or it will be 0. lambda = (double)1/(N_cluster_ScaledBack+1);
Fixed Bug in LogLikelihoodCluster and LogLikelihoodNonCluster: for synonymous, recurrent is not considered; or the replacement recurrent M will be calculated, and create a bug.
*/

#include "cPRFCluster.h"
#include <string>
#include <list>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <thread>

using namespace std;
#ifdef __linux__
unsigned int NUM_CPU = get_nprocs_conf();
#else
unsigned int NUM_CPU = std::thread::hardware_concurrency();
#endif

/***************************************************
* Function: Initialization of parameters
***************************************************/
//Declare and initiate variables for the function cPRFCluster
cPRFCluster::cPRFCluster() {
  flag_N_pol=0; //polymorphism
  flag_N_div=0; //divergence

  confidence_interval=0.95;
  quantile_for_CI=(1.0-confidence_interval)/2.0; //0.025

  flag_found_pr=0; //polymorphism replacement
  flag_found_dr=0; //divergence replacement
  flag_found_ps=0; //polymorphism sysnonymous
  flag_found_ds=0; //divergence synonymous

  //Initiate parameters, require inputs from the user
  pol_cons_seqfile = ""; // required user input with the option '-pc file_name', the polymorphism file name with sequences in format of R, S, *.
  div_cons_seqfile = ""; // required user input with the option '-dc file_name', the divergence file name with sequences in format of R, S, *.

  SilentRate = 0.0;
  SilentRate_flag=0;

  polymorphism_num_s =""; // required user input with the option '-pn polymorphism_seq_number'.
  tumor_num_s =""; // required user input with the option '-dn divergence_seq_number'.

  polymorphism_num =-299;
  tumor_num =-299;
  output_format_num=0;
  genetic_code = 1;
  criterion_type = 0;
  Do_Synonymous_Cluster=0; // By default, option for clustering Synonymous sites in the polymorphism and divergence sequence is off. Users need to use '-s 1' to turn it on.
  MS_only=0;
  ci_ma=0; // confidence interval for model average
  Div_time=0.0; // Species divergent time; can be initiate by the user using prior species divergent time, or it will be estimated based on the given sequence, and it may be biased by the gene.
  Do_r_estimate=1; // variable for the estimated selection coefficient r.
  Do_ci_r=1; // confidence intervals for selection coefficient
  Do_ci_r_exact=0; //exact algorithm estimated r confidence interval
  //ZMZ 04/28/2016  added option - regional gamma only
  regional_gamma_only=0; //regional only gamma option when regional_gamma_only=1, default weighted gamma
  Nuc_replace=1;
  NI_estimate=0; // By default, Nuetrality Index will not be estimated. The user can turn it on with '-NI 1'.
  TotalRecurrentCount=0;
  TotalRecurrentSite=0;
  TotalReplacementSite=0;
  TotalSilentSite=0;
}

/***************************************************
* Function: remove elements in the vector
***************************************************/
cPRFCluster::~cPRFCluster() {

  div_cons_seq.clear();
  div_cons_seqname.clear();

  vec_r_c.clear();
  vec_r_c_r.clear();
  vec_rModels_c.clear();
  vec_lower_r_c.clear();
  vec_upper_r_c.clear();
  
  vec_lower_r_c_r.clear();
  vec_upper_r_c_r.clear();

  vec_lower_rate_ds.clear();
  vec_lower_rate_dr.clear();
  
  vec_upper_rate_ds.clear();
  vec_upper_rate_dr.clear();

  vec_MA_rate_ds.clear();
  vec_MA_rate_dr.clear();
  
  vec_MS_rate_ds.clear();
  vec_MS_rate_dr.clear();

  vec_SelectedModels.clear();
  vec_SelectedModels_ds.clear();
  vec_SelectedModels_dr.clear();
  
  vec_AllModels.clear();
  vec_AllModels_ds.clear();
  vec_AllModels_dr.clear();

}

/***************************************************
* Function: initialization; change the size of the vectors based on gene length and empty the vector
***************************************************/
int cPRFCluster::init(long N){
  vec_r_c.resize(N,0.0);
  vec_r_c_r.resize(N,0.0);
  //vec_rModels_c.resize(N,0.0);
  vec_lower_r_c.resize(N,0.0);
  vec_upper_r_c.resize(N,0.0);

  vec_lower_r_c_r.resize(N,0.0);
  vec_upper_r_c_r.resize(N,0.0);

  vec_MS_rate_ds.resize(N,0.0);
  vec_MS_rate_dr.resize(N,0.0);

  vec_MA_rate_ds.resize(N,0.0);
  vec_MA_rate_dr.resize(N,0.0);

  vec_lower_rate_ds.resize(N,0.0);
  vec_lower_rate_dr.resize(N,0.0);

  vec_upper_rate_ds.resize(N,0.0);
  vec_upper_rate_dr.resize(N,0.0);

  vec_SelectedModels_ds.clear();
  vec_SelectedModels_dr.clear();

  vec_AllModels_ds.clear();
  vec_AllModels_dr.clear();

  vec_SelectedModels.clear();
  vec_AllModels.clear();
  vec_MA_rate.resize(N,0.0);
  vec_lower_rate.resize(N,0.0);
  vec_upper_rate.resize(N,0.0);
  return 1;
}
/***************************************************
* Function: Read the input files and execute the main function RunML, and screen output;
* Input Parameter: required divergence file names, and the number of sequences in divergence
***************************************************/
int cPRFCluster::Run(int argc, const char*argv[]) {	
  int i, flag=1;
  srand(1234); // to fix random number generator, to make sure the program is repeatable; need to remove for the final version of the program
  try {		
	    //Write in the output file with program name, Version, LastUpdate, and Reference.
	    cout<<endl<<NAME<<", Version: "<<VERSION<<" [Last Update: "<<LASTUPDATE<<"]"<<endl;
	    cout<<"Reference: "<<REFERENCE<<endl<<endl;
	    static time_t time_start = time(NULL); // Record the start time
	  //Parse input parameters
    if (parseParameter(argc, argv)!=1) throw "Error in parsing parameters!";
    
    //Check for presence of the input file and whether the number of tumors is specified
    if(div_cons_seqfile=="" || tumor_num_s=="") throw "Failed to specify the input file and/or the number of tumors! Use -H to find out.";
    //convert the format for the sequence number from string to int
    tumor_num=CONVERT<int>(tumor_num_s);
    cout<<"Read cancer divergence consensus input file: "<<div_cons_seqfile<<endl<<endl;
    if (readFasta(div_cons_seqfile, div_cons_seqname, div_cons_seq)!=1) throw "Error in reading divergent sequence.";
    //Print the divergence consensus sequence and the gene name
    cout<<endl<<"Divergence Consensus Sequence:"<<endl<<">"<<div_cons_seqname[0].c_str()<<endl;
    cout<<div_cons_seq[0].c_str()<<endl<<endl<<endl;
    //Check sequence length for divergence
   if(div_cons_seq[0].size()%3!=0) cout<< "Warning: the length of the specified divergence sequence cannot be divided by 3 (codon size)."<<endl;

   //Get the recurrent position and count in Recurrents by taking recurrent.txt as an input as the seq file, flag
   cout<<"Recurrent positions and counts: "<<endl;
   GetRecurrentList(recurrent_file); // JT: rename as GetRecurrentList - ZMZ fixed
   cout<<"Number of recurrent positions: "<<Recurrents.size()<<endl;
   long jj;
   for (jj=0;jj<Recurrents.size();jj++){
	   cout<<"Recurrent Site:\t"<<Recurrents[jj].sites<<"\tCount: "<<Recurrents[jj].counts<<endl;
	   int count=Recurrents[jj].counts;
	   TotalRecurrentCount+=count;
	   TotalRecurrentSite+=1;
   }
   //Calculate the gamma and 95% CI gamma for each recurrent site, using 2r/(1-e^(-2r))=RecurrentNumber/(ReplacementRate*TumorNumber)
   //Get the lookup table read for CI for all different k
   //ZMZ 06/16/2016 update LookupTable_CSIMAC_CI_Recurrent_v10.dat from v9, now with 1000 recur counts
   string input_lookup_file="LookupTable_CSIMAC_CI_Recurrent_v10.dat";
   LambdaCIs.clear();
   LambdaCILookupTable(input_lookup_file);
   cout<<"Maximum recurrent count currently permitted by CSIMAC: "<<LambdaCIs.size()<<endl;
   if (LambdaCIs.size()==0) 
   {   
    cout<<"Error: The LambdaCILookupTable "<<input_lookup_file<<" is empty...\n"; // JT: fix this with a stringconcat or something - ZMZ fixed
    throw;   
   }

   //**** Main Step: Run the main function Maximum Likelihood for cancer divergence sequences, to get site specific gamma and 95% CI gamma
   RunML(div_cons_seq);

    //Display on screen after finish running the program and print out the time used.
    cout<<endl<<"Mission accomplished. (Time elapsed: ";
    time_t t = time(NULL)-time_start;
    int h=t/3600, m=(t%3600)/60, s=t-(t/60)*60;
    if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
    else   cout<<m<<":"<<s<<")"<<endl;
  }	
  catch (const char* e) {
    cout<<e<<endl;		
    flag = 0;
  }
  catch (...) {		
    flag = 0;
  }
  return flag;
}
/***************************************************
* Function: Print out the cluster information and gamma values for replacement divergence by default; and information for synonymous divergence can be printed out if the user assign '-s 1'.
***************************************************/
int cPRFCluster::output(long N){
   cout<<endl<<"//Results based on model selection: "<<endl;
   cout<<"Cancer divergence synonymous mutation rate (ucs): "<<ucs<<endl;
   cout<<"Cancer divergence replacement mutation rate (ucr): "<<ucr<<endl;
  // Print out cluster information for synonymous sites in the Divergence sequence if the user assigns '-s 1'.
  if(Do_Synonymous_Cluster==1){
    cout<<endl<<"Clusters of synonymous divergence:"<<endl;
    if(vec_SelectedModels_ds.size()==0){
      cout<<"Note: Synonymous divergence (DS) = 1 or 0. There are not enough synonymous divergent sites for clustering!"<<endl<<endl;
    }else if(vec_SelectedModels_ds.size()==1 && vec_SelectedModels_ds[0].pos_start==vec_SelectedModels_ds[0].cs && vec_SelectedModels_ds[0].pos_end==vec_SelectedModels_ds[0].ce){
      cout<<"Note: CSIMAC identified no clustering of synonymous sites in this gene."<<endl<<endl;
    }else{
      // JT: to debug vec_SelectedModels_dr and vec_SelectedModels_ds output models. - ZMZ fixed.
      for(long i=0; i<vec_SelectedModels_ds.size(); i++){
    	  if (output_format_num==1)
    	{
    		cout<<"\nNucleotide\nStart_Position = "<<(vec_SelectedModels_ds[i].pos_start*Scale+1)<<"\tEnd_Position = "<<(vec_SelectedModels_ds[i].pos_end*Scale+1);
    		cout<<"\tCluster_Start_Position = "<<(vec_SelectedModels_ds[i].cs*Scale+1)<<"\tCluster_End_Position = "<<(vec_SelectedModels_ds[i].ce*Scale+1);
    	}
    	else if (output_format_num==0)
    	{
    		cout<<"\nAmino Acid\nStart_Position = "<<(vec_SelectedModels_ds[i].pos_start*Scale/3+1)<<"\tEnd_Position = "<<(vec_SelectedModels_ds[i].pos_end*Scale/3+1);
    		cout<<"\tCluster_Start_Position= "<<(vec_SelectedModels_ds[i].cs*Scale/3+1)<<"\tCluster_End_Position= "<<(vec_SelectedModels_ds[i].ce*Scale/3+1);
    	}
	cout<<endl;	
	cout<<"InL_0= "<<vec_SelectedModels_ds[i].InL0<<"\tInL= "<<vec_SelectedModels_ds[i].InL;
	cout<<"\tAIC_0= "<<vec_SelectedModels_ds[i].AIC0<<"\tAIC= "<<vec_SelectedModels_ds[i].AIC;
	cout<<"\tAICc_0= "<<vec_SelectedModels_ds[i].AICc0<<"\tAICc= "<<vec_SelectedModels_ds[i].AICc;
	cout<<"\tBIC_0= "<<vec_SelectedModels_ds[i].BIC0<<"\tBIC= "<<vec_SelectedModels_ds[i].BIC;
	cout<<endl;	
	cout<<"P0_DivergenceSynonymous= "<<vec_SelectedModels_ds[i].p0<<"\tPc_DivergenceSynonymous= "<<vec_SelectedModels_ds[i].pc;
	cout<<endl<<endl;
      }
    }
  }

  // Print out cluster information for Replacement sites in the Divergence sequence by default.
  cout<<endl<<"Clusters of replacement divergence:"<<endl;
  if(vec_SelectedModels_dr.size()==0){
    cout<<"Note: Replacement divergence (DR) = 0. There is no replacement divergent sites for clustering!"<<endl<<endl;
  }else if(vec_SelectedModels_dr.size()==1 && vec_SelectedModels_dr[0].pos_start==vec_SelectedModels_dr[0].cs && vec_SelectedModels_dr[0].pos_end==vec_SelectedModels_dr[0].ce){
    cout<<"Note: There is no cluster of replacement sites in this gene."<<endl<<endl;
  }else{  
    // JT: to debug vec_SelectedModels_dr and vec_SelectedModels_ds output models. - ZMZ fixed
    for(long i=0; i<vec_SelectedModels_dr.size(); i++){ 

  	  if (output_format_num==1) //0: amino acid output || 1: nucleotide output, default=0
  	{
  		cout<<"\nNucleotide\nStart_Position = "<<(vec_SelectedModels_dr[i].pos_start*Scale+1)<<"\tEnd_Position = "<<(vec_SelectedModels_dr[i].pos_end*Scale+1);
  		cout<<"\tCluster_Start_Position= "<<(vec_SelectedModels_dr[i].cs*Scale+1)<<"\tCluster_End_Position= "<<(vec_SelectedModels_dr[i].ce*Scale+1);
  	}
  	else if (output_format_num==0) //0: amino acid output || 1: nucleotide output, default=0
  	{
  		cout<<"\nAmino Acid\nStart_Position = "<<(vec_SelectedModels_dr[i].pos_start*Scale/3+1)<<"\tEnd_Position = "<<(vec_SelectedModels_dr[i].pos_end*Scale/3+1);
  		cout<<"\tCluster_Start_Position= "<<(vec_SelectedModels_dr[i].cs*Scale/3+1)<<"\tCluster_End_Position= "<<(vec_SelectedModels_dr[i].ce*Scale/3+1);
  	}
      cout<<endl;
      cout<<"InL_0= "<<vec_SelectedModels_dr[i].InL0<<"\tInL= "<<vec_SelectedModels_dr[i].InL;
      cout<<"\tAIC_0= "<<vec_SelectedModels_dr[i].AIC0<<"\tAIC= "<<vec_SelectedModels_dr[i].AIC;
      cout<<"\tAICc_0= "<<vec_SelectedModels_dr[i].AICc0<<"\tAICc= "<<vec_SelectedModels_dr[i].AICc;
      cout<<"\tBIC_0= "<<vec_SelectedModels_dr[i].BIC0<<"\tBIC= "<<vec_SelectedModels_dr[i].BIC;
      cout<<endl;
      cout<<"P0_DivergenceReplacement= "<<vec_SelectedModels_dr[i].p0<<"\tPc_DivergenceReplacement= "<<vec_SelectedModels_dr[i].pc;
      cout<<endl<<endl;
    }
  }

  if (N*Scale%3!=0) { cout<< "Warning: the length of divergence sequence cannot be divided by 3 (codon size)."<<endl;}

  //Print the results of Model Averaging for each position
  if (MS_only==0) {
    cout<<endl<<"//Results based on model averaging: "<<endl;
    cout.setf(ios::left);
    int width=15;

    //Output the title
    cout.width(width); cout<<"Position\t";
    if(Do_Synonymous_Cluster==1){ // JT: could we change to Do_Syn_cluster or Do_SynSite_Cluster or Do_Synonymous_Cluster? - ZMZ fixed
        cout.width(width); cout<<"MS_DivRep\t";
        cout.width(width); cout<<"MA_DivRep\t";
    	cout.width(width); cout<<"MS_DivSys\t";
        cout.width(width); cout<<"MA_DivSys\t";
      if(ci_ma==1){
          cout.width(width); cout<<"Lower_CI_DivRep\t";
          cout.width(width); cout<<"Upper_CI_DivRep\t";
    	  cout.width(width); cout<<"Lower_CI_DivSys\t";
	      cout.width(width); cout<<"Upper_CI_DivSys\t";
      }
    }

    if(Do_r_estimate==1){ // JT: change to Do_r_estimate - ZMZ fixed
      cout.width(width); cout<<"Gamma_Cancer";
      if(Do_ci_r==1){ // JT: change do Do_ci_r - ZMZ fixed
		cout.width(width); cout<<"\tLower_CI_Gamma_c\t";
		cout.width(width); cout<<"Upper_CI_Gamma_c";
		cout.width(width); cout<<"\tMutSymbol_Q2_D-2_R1_S-1_*0\t"; //Q stands for Recurrent using 2; D stands for Damaging using -2; R stands for Replacement using 1; S stands for Silent using -1; * stands for conserved using 0
		cout.width(width); cout<<"MutationStatus"; // MutationStatus: Q (recurrent), D (Damaging), R (replacement), S (silent) or *(no mutation); 
      }
    }
    cout<<endl;//End of the title line of the output values.
    //Output the data in the format of nucleotide or amino acid sequence, 0: amino acid output || 1: nucleotide output, default=0
    if (output_format_num==1 || (output_format_num==0 && Scale>1)) //(output_format_num==0 and Scale==3) and (output_format_num==1 and Scale==1) are exactly the same output position
    {
		if (output_format_num==0 and Scale>1) {Scale=Scale/3;}//for amino acid, scale it down 3 times to make it comparable to nucleotide Scale
		for(long i=0; i<N; i++) {//each position after scaled
			for (long j=0;j<Scale;j++){ // JT: fix this Scale/3 = 1/3 - ZMZ fixed adjusting amino acid scales down 3 times to be comparable to nucleotide and expanding each position i to i*Scale+j (0<j<Scale).
				cout.width(width);cout<<Scale*i+j+1<<"\t"; //convert to the real amino acid or nucleotide position of the gene, if Scale==1, i; if Scale=6, i=10, new position would be 60-65. 
				if(Do_Synonymous_Cluster==1){
					cout.width(width);cout<<vec_MS_rate_dr[i]<<"\t";
					cout.width(width);cout<<vec_MA_rate_dr[i]<<"\t";
					cout.width(width);cout<<vec_MS_rate_ds[i]<<"\t";
					cout.width(width);cout<<vec_MA_rate_ds[i]<<"\t";
					if(ci_ma==1){
						cout.width(width);cout<<vec_lower_rate_dr[i]<<"\t";
						cout.width(width);cout<<vec_upper_rate_dr[i]<<"\t";
						cout.width(width);cout<<vec_lower_rate_ds[i]<<"\t";
						cout.width(width);cout<<vec_upper_rate_ds[i]<<"\t";
					}
				  }
				//Do_r_estimate: choose the option of outputting gamma, and 95% confidence intervals of gamma
				if(Do_r_estimate==1){
					cout.width(width); //cancer gamma
					//gamma for site i
					if(vec_r_c[i]==299){
					  cout.width(width);cout<<"INF";
					}else if (vec_r_c[i]==-299){
					  cout.width(width);cout<<"N-INF";
					}else if (vec_r_c[i]==-199){
					  cout.width(width);cout<<"NULL";
					}else{
					  cout.width(width);cout<<vec_r_c[i];
					}
					  //Lower confidence interval of gamma for site i
					  if(vec_lower_r_c[i]==299){
						cout<<"\t";cout.width(width);cout<<"INF"<<"\t";
					  }else if(vec_lower_r_c[i]==-299){
						cout<<"\t";cout.width(width);cout<<"N-INF"<<"\t";
					  }else if(vec_lower_r_c[i]==-199){
						cout<<"\t";cout.width(width);cout<<"NULL"<<"\t";
					  }else{
						cout<<"\t";cout.width(width);cout<<vec_lower_r_c[i]<<"\t";
					  }
					  //Upper confidence interval of gamma for site i
					  if(vec_upper_r_c[i]==299){
						cout.width(width);cout<<"INF";
					  }else if(vec_upper_r_c[i]==-299){
						cout.width(width);cout<<"N-INF";
					  }else if(vec_upper_r_c[i]==-199){
						cout.width(width);cout<<"NULL";
					  }else{
						cout.width(width);cout<<vec_upper_r_c[i];
						}
				  }//end of if(Do_r_estimate==1)		  
				//Site labels using number 0, -1, 1, 2, 3
				  if (div_codon_consensus[i]=='*') { cout.width(width);cout<<"\t"<<0; }
				  else if (div_codon_consensus[i]=='S') { cout.width(width);cout<<"\t"<<-1; } //Re-ordered
				  else if (div_codon_consensus[i]=='R') { cout.width(width);cout<<"\t"<<1; }
				  else if (div_codon_consensus[i]=='Q') { cout.width(width);cout<<"\t"<<2; } //recurrent sites
				  else if (div_codon_consensus[i]=='D') { cout.width(width);cout<<"\t"<<-2; } //Added Damaging mutation records
				  else { throw 1;}
				  cout<<"\t"<<div_codon_consensus[i];
				  cout<<endl;
			}//end of the inside for loop
		  }//end of the outside for loop
	  cout<<endl;
    }//end of if (output_format_num==1 or (output_format_num==0 and Scale!=1)) 
  
  //Output the data in the format of amino acids without scaling
  if (output_format_num==0 and Scale==1) // 0: amino acid output || 1: nucleotide output, default=0
  {
	  for(long i=0; i<N/3; i++) { 
		   cout.width(width);cout<<i+1<<"\t";
		  if(Do_Synonymous_Cluster==1){
			    cout.width(width);cout<<(vec_MS_rate_dr[i*3]+vec_MS_rate_dr[i*3+1]+vec_MS_rate_dr[i*3+2])/3<<"\t";
			    cout.width(width);cout<<(vec_MA_rate_dr[i*3]+vec_MA_rate_dr[i*3+1]+vec_MA_rate_dr[i*3+2])/3<<"\t";
			    cout.width(width);cout<<(vec_MS_rate_ds[i*3]+vec_MS_rate_ds[i*3+1]+vec_MS_rate_ds[i*3+2])/3<<"\t";
				cout.width(width);cout<<(vec_MA_rate_ds[i*3]+vec_MA_rate_ds[i*3+1]+vec_MA_rate_ds[i*3+2])/3<<"\t";
				if(ci_ma==1){
					cout.width(width);cout<<(vec_lower_rate_dr[i*3]+vec_lower_rate_dr[i*3+1]+vec_lower_rate_dr[i*3+2])/3<<"\t";
					cout.width(width);cout<<(vec_upper_rate_dr[i*3]+vec_upper_rate_dr[i*3+1]+vec_upper_rate_dr[i*3+2])/3<<"\t";
					cout.width(width);cout<<(vec_lower_rate_ds[i*3]+vec_lower_rate_ds[i*3+1]+vec_lower_rate_ds[i*3+2])/3<<"\t";
					cout.width(width);cout<<(vec_upper_rate_ds[i*3]+vec_upper_rate_ds[i*3+1]+vec_upper_rate_ds[i*3+2])/3<<"\t";
				}
			}//end of if(Do_Synonymous_Cluster==1)
			//Do_r_estimate: choose the option of outputting gamma, and 95% confidence intervals of gamma
			if(Do_r_estimate==1){
				cout.width(width);
				//gamma for site i
				if(vec_r_c[i*3]==299 or vec_r_c[i*3+1]==299 or vec_r_c[i*3+2]==299){
				  cout.width(width);cout<<"INF";
				}else if (vec_r_c[i*3]==-299 or vec_r_c[i*3+1]==-299 or vec_r_c[i*3+2]==-299){
				  cout.width(width);cout<<"N-INF";
				}else if (vec_r_c[i*3]==-199 || vec_r_c[i*3+1]==-199 || vec_r_c[i*3+2]==-199){
				  cout.width(width);cout<<"NULL";
				}else{
				  cout.width(width);cout<<(vec_r_c[i*3]+vec_r_c[i*3+1]+vec_r_c[i*3+2])/3;
				}
				 //Lower confidence interval of gamma for site i
				  if(vec_lower_r_c[i*3]==299 or vec_lower_r_c[i*3+1]==299 or vec_lower_r_c[i*3+2]==299){
					cout<<"\t";cout.width(width);cout<<"INF"<<"\t";
				  }else if(vec_lower_r_c[i*3]==-299 or vec_lower_r_c[i*3+1]==-299 or vec_lower_r_c[i*3+2]==-299){
					cout<<"\t";cout.width(width);cout<<"N-INF"<<"\t";
				  }else if(vec_lower_r_c[i*3]==-199 || vec_lower_r_c[i*3+1]==-199 || vec_lower_r_c[i*3+2]==-199){
					cout<<"\t";cout.width(width);cout<<"NULL"<<"\t";
				  }else{
					cout<<"\t";cout.width(width);cout<<(vec_lower_r_c[i*3]+vec_lower_r_c[i*3+1]+vec_lower_r_c[i*3+2])/3<<"\t";
				  }
				  //Upper confidence interval of gamma for site i
				  if(vec_upper_r_c[i*3]==299 or vec_upper_r_c[i*3+1]==299 or vec_upper_r_c[i*3+2]==299){
					cout.width(width);cout<<"INF";
				  }else if(vec_upper_r_c[i*3]==-299 or vec_upper_r_c[i*3+1]==-299 or vec_upper_r_c[i*3+2]==-299){
					cout.width(width);cout<<"N-INF";
				  }else if(vec_upper_r_c[i*3]==-199 || vec_upper_r_c[i*3+1]==-199 || vec_upper_r_c[i*3+2]==-199){
					cout.width(width);cout<<"NULL";
				  }else{
					cout.width(width);cout<<(vec_upper_r_c[i*3]+vec_upper_r_c[i*3+1]+vec_upper_r_c[i*3+2])/3;
					}
			  }	//end of if(Do_r_estimate==1)		  
			//Site labels using number 0, -1, 1, 2, 3
			if (div_codon_consensus[i*3]=='*' and div_codon_consensus[i*3+1]=='*' and div_codon_consensus[i*3+2]=='*') { cout.width(width);cout<<"\t"<<0; }
			else if (div_codon_consensus[i*3]=='S' or div_codon_consensus[i*3+1]=='S' or div_codon_consensus[i*3+2]=='S') { cout.width(width);cout<<"\t"<<-1; } // Reordered; S first, can be overwritten by laters
			else if (div_codon_consensus[i*3]=='R' or div_codon_consensus[i*3+1]=='R' or div_codon_consensus[i*3+2]=='R') { cout.width(width);cout<<"\t"<<1; }
			else if (div_codon_consensus[i*3]=='Q' or div_codon_consensus[i*3+1]=='Q' or div_codon_consensus[i*3+2]=='Q') { cout.width(width);cout<<"\t"<<2; }
			else if (div_codon_consensus[i*3]=='D' or div_codon_consensus[i*3+1]=='D' or div_codon_consensus[i*3+2]=='D') { cout.width(width);cout<<"\t"<<-2; } // Added Damaging mutation records
			else { throw 1;}
			//Site labels using characters *, S, R, Q, D
			if (div_codon_consensus[i*3]=='*' and div_codon_consensus[i*3+1]=='*' and div_codon_consensus[i*3+2]=='*') { cout<<"\t*"; }
			else if (div_codon_consensus[i*3]=='S' or div_codon_consensus[i*3+1]=='S' or div_codon_consensus[i*3+2]=='S') {	cout<<"\tS";  }
			else if (div_codon_consensus[i*3]=='R' or div_codon_consensus[i*3+1]=='R' or div_codon_consensus[i*3+2]=='R') { cout<<"\tR";  }		    
			else if (div_codon_consensus[i*3]=='Q' or div_codon_consensus[i*3+1]=='Q' or div_codon_consensus[i*3+2]=='Q') { cout<<"\tQ";  }
			else if (div_codon_consensus[i*3]=='D' or div_codon_consensus[i*3+1]=='D' or div_codon_consensus[i*3+2]=='D') {	cout<<"\tD";  }
			else { throw 1;}
			cout<<endl;//end of site i, go to next site in the for loop
		  }//end of the for loop, for(long i=0; i<N/3; i++)
  }//end of if (output_format_num==0 and Scale==1) 
}// end of if (MS_only==0)
else if (MS_only==1){
    cout<<endl<<"*************"<<endl;
    cout<<"Warning:"<<endl<<"Check the parameter -m. If only model selection is used to find the probability of the site being a variant is performed, there is no estimate of the selection coefficient (gamma) and its confidence intervals. Please check tutorial for more details!"<<endl;
    cout<<"*************"<<endl;
    return 1;
}
  cout<<endl<<"Abbreviation: MS=Model Selection; MA=Model Averaging; CI=Confidence Interval; ds=Divergence Synonymous; dr=Divergence Replacement; Gamma=N*s (Gamma: scaled selection coefficient (selection intensity); s: selection coefficient); gamma >1 Negative selection, <1 Positive selection); INF=Infinite; N-INF=Negative Infinite; NULL=Not enough information for this site"<<endl;
  cout<<"Abbreviation: MutSymbol_Q2_D-2_R1_S-1_*0: Q stands for Recurrent using 2; D stands for Damaging using -2; R stands for Replacement using 1; S stands for Silent using -1; * stands for conserved using 0."<<endl;
  cout<<endl<<"#End of clustering"<<endl<<endl;
  return 1;
}

/***************************************************
* Function: Main function for Clustering by Maximum likelihood for synonymous and replacement sites in the polymorphism and divergence sequences.
* Input Parameter: polymorphism sequence, divergence sequence
* Output:
* Return Value:
***************************************************/
int cPRFCluster::RunML(vector<string> div_cons_seq) {
   //use only one format of input for the final version
   div_codon_consensus = div_cons_seq[0]; // Get divergence sequence with synonymous (S) and replacement (R) sites labeled 
   long N=div_codon_consensus.length(); //polymorphism and divergence sequence length, the two are equal.
   init(N); //initialization subfunction: change the size of the vectors based on gene length and empty the vector
   GeneLength=N;

  double ds=0.0;// initializing synonymous divergence at zero
  double dr=0.0;// initializing replacement divergence at zero
 
  // Count synonymous divergent sites by finding symbol 'S' in the divergent sequence from start position 0 to the end position N-1
  ds=getDifference(div_codon_consensus,0,N-1,'S');
  // Count replacement divergent sites by finding symbol 'R' in the divergent sequence from start position 0 to the end position N-1
  dr=getDifference(div_codon_consensus,0,N-1,'R');
  TotalReplacementSite=dr;
  TotalSilentSite=ds;

	//Print the # of divergent synonymous and replacement sites labeled with 'S' and 'R'
  cout<<"Synonymous divergent count (DS): "<<ds<<endl;
  cout<<"Replacement divergent count (DR): "<<dr<<endl;
  cout<<"Number of tumor samples: "<<tumor_num<<endl;
  long N_ScaledBack=N*Scale; //Scale the gene length back 
  cout<<"Gene length before scaling: "<<GeneLength<<"\tAfter scaling: "<<N_ScaledBack<<" bp"<<endl;

//Get the synonymous mutation rate (ucs) based on the number of silent mutations in the gene or from user-input (MutSigCV SilentRate). The latter is preferred.
  ucs=CancerSynonymousRate(tumor_num, N);

   //double ratio_NS=CalculateNS(ref_seq_gene); // calculate the ratio of replacement and synonymous N/S
   double ratio_NS= 0.345291479; //the ratio is (0.77/2.23), assuming that the first and third position will be nonsynonymous 5% and 72% of the time. See Nei and Gojobori's paper 1986
   cout<<"The ratio of replacement mutations to synonymous mutations is "<<ratio_NS<<endl;
  ucr=ReplacementRate(ratio_NS, ucs); // Get replacement rate for cancer divergence (ucr)
  cout<<"Cancer replacement divergence mutation rate: "<<ucr<<endl;

	// Fixed bug due to the new OS X 10.9 system [error: variable length array of non-POD element type 'struct SiteModels']. Solution: use a very large number instead of the parameter N for the gene length, for keeping all models for each gene site, to make sure the number is larger than the gene length.
  // struct SiteModels sm_pol[N];
  struct SiteModels sm_div[10000];
  if (N>10000) { cout<<"The length of the gene exceeds 10000. Scale it down or revise the SiteModels upper-boundary array size!"<<endl; throw 1;}

  //cout<<" Estimated time: "<<(12*N/1000+8)<<" minutes for this gene with 20000 models."<<endl;
  cout<<"CSIMAC is conducting the two time-consuming steps, ClusterSubSeq and Do_ci_r_stochastic."<<endl;
  //cout<<" ClusterSubSeq depends on the model number, rate: ~4 minutes for each additional 10000 models, ~6 hours for 1 million models"<<endl;
  //cout<<" Do_ci_r_stochastic is determined by gene length, with the speed rate 12 minutes per 1000 sites."<<endl;

  // Clustering and 95% CI for the clusters, this will be calculated for Silent divergence if (Do_Synonymous_Cluster==1)
  if(Do_Synonymous_Cluster==1){
    //Initialize for DS
    vec_SelectedModels.clear();
    vec_MS_rate.clear();
    vec_MA_rate.clear();
    vec_lower_rate.clear();
    vec_upper_rate.clear();

    vec_MS_rate.resize(N,0.0);
    vec_MA_rate.resize(N,0.0);
    vec_lower_rate.resize(N,0.0);
    vec_upper_rate.resize(N,0.0);

    vec_MS_rate_ds.clear();
    vec_MA_rate_ds.clear();
    vec_MS_rate_ds.resize(N,0.0);
    vec_MA_rate_ds.resize(N,0.0);

	// Major Step: Find cluster and calculate probability using multiple models for synonymous divergence
    cout<<endl<<"Starting the clustering of synonymous variants"<<endl;
    ClusterSubSeq(0, N-1,'S',sm_div); // Major subroutine for clustering of synonymous variants
    cout<<"Finished clustering of synonymous variants."<<endl;
    vec_SelectedModels_ds=vec_SelectedModels;
    vec_MS_rate_ds=vec_MS_rate;
    vec_MA_rate_ds=vec_MA_rate;
    if(MS_only==0 && ci_ma==1 && ds >0){
      CI_MA(sm_div,N); // Major subroutine for 95% CI for probability of synonymous variant for each site
    }
    vec_lower_rate_ds=vec_lower_rate;
    vec_upper_rate_ds=vec_upper_rate;
  }

  //Empty vectors and re-size for replacement divergence (DR)
  vec_SelectedModels.clear();
  vec_SelectedModels_dr.clear();
  vec_MS_rate.clear();
  vec_MA_rate.clear();
  vec_lower_rate.clear();
  vec_upper_rate.clear();

  vec_MS_rate.resize(N,0.0);
  vec_MA_rate.resize(N,0.0);
  vec_lower_rate.resize(N,0.0);
  vec_upper_rate.resize(N,0.0);

  vec_MS_rate_dr.clear();
  vec_MA_rate_dr.clear();
  vec_MS_rate_dr.resize(N,0.0);
  vec_MA_rate_dr.resize(N,0.0);

	//Major Step: Find cluster and calculate probability using multiple models for replacement divergence
  cout<<endl<<"Starting the clustering of replacement variants"<<endl;
  time_t time_start1 = time(NULL); // Record the start time
  ClusterSubSeq(0, N-1,'R',sm_div); //find clusters in the sequence
  cout<<"Finished clustering of replacement variants."<<endl;
  cout<<"Time elapsed during clustering of replacement variants: ";
  time_t t2 = time(NULL)-time_start1;
  int h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
  else   cout<<m<<":"<<s<<")"<<endl;

  vec_SelectedModels_dr=vec_SelectedModels; // all selected models for the clustering
  vec_MS_rate_dr=vec_MS_rate; // selected model rate for the clustering rate
  vec_MA_rate_dr=vec_MA_rate; // model averaged rate for the clustering rate

  // Calculate the 95% CI for the clustering probability
  if(MS_only==0 && ci_ma==1 && dr>0){
	  cout<<endl<<"Starting model-averaging analysis to yield CIs for replacement divergence"<<endl;
	  time_t time_start1 = time(NULL); // Record the start time
      CI_MA(sm_div,N);
	  cout<<endl<<"Finished model-averaging analysis to yield CIs for replacement divergence"<<endl;
	  cout<<"Time elapsed during model-averaging analysis to yield CIs: ";
	  time_t t2 = time(NULL)-time_start1;
	  int h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
	  if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
	  else   cout<<m<<":"<<s<<")"<<endl;
  }
  vec_lower_rate_dr=vec_lower_rate; // lower 95% CI rate for the clustering probability
  vec_upper_rate_dr=vec_upper_rate; // upper 95% CI rate for the clustering probability

 //Estimate gamma for Replacement sites
  if (MS_only==0 && Do_r_estimate==1) {
    //If dr == 0, r can't be estimated.
    if(dr>0){
      cout<<endl<<"Starting estimation of gamma for human cancer replacement divergence."<<endl;
      time_t time_start1 = time(NULL); // Record the start time
      rc_SitePRF(tumor_num,ucr, N); //estimate gamma for cancer divergence
      cout<<"Estimation of gamma for human cancer replacement divergence completed."<<endl;
      cout<<"Time elapsed during estimation of gamma for human cancer replacement divergence: ";
      time_t t2 = time(NULL)-time_start1;
      int h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
      if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
      else   cout<<m<<":"<<s<<")"<<endl;
    }else{
      cout<<endl<<"*************"<<endl;
      cout<<"There are no replacement sites to estimate gamma!"<<endl;
      cout<<"*************"<<endl;
    }

    // Calculate confidence intervals for r using the exact or stochastic algorithm.
    if(dr>0 && Do_ci_r==1){
      if(Do_ci_r_exact==1){ // Use exact algorithm to calculate confidence intervals for r
        cout<<endl<<"Starting the estimation of CIs for Gamma by the exact algorithm."<<endl;
        time_t time_start1 = time(NULL); // Record the start time
        CIr_exact(sm_div,N);
        cout<<endl<<"Estimation of CIs for Gamma by exact algorithm complete."<<endl;
        cout<<"Time elapsed during estimation of CIs for Gamma by exact algorithm: ";
          time_t t2 = time(NULL)-time_start1;
          int h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
          if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
          else   cout<<m<<":"<<s<<")"<<endl;
      }else{ // Use stochastic algorithm to calculate r Confidence Intervals
        cout<<endl<<"Starting the estimation of CIs for Gamma by stochastic algorithm."<<endl;
        time_t time_start1 = time(NULL); // Record the start time
        CIr_stochastic(sm_div,N);
        cout<<endl<<"Estimation of CIs for Gamma by stochastic algorithm complete."<<endl;
        cout<<"CIr_stochastic (Time elapsed: ";
          time_t t2 = time(NULL)-time_start1;
          int h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
          if(h)  cout<<h<<":"<<m<<":"<<s<<")"<<endl;
          else   cout<<m<<":"<<s<<")"<<endl;
      }
    }
    //for cases replacement sites  and divergence sequences are too few.
    else if(Do_ci_r==1){
      cout<<endl<<"*************"<<endl;
      cout<<"There are no replacement sites to estimate CIs for gamma!"<<endl;
      cout<<"*************"<<endl;
    }    
  }//end of if (MS_only==0 && Do_r_estimate==1)
       
  if(MS_only==1){
    cout<<endl<<"*************"<<endl;
    cout<<"Warning:"<<endl<<"Check the parameter -m. If only model selection is used to find the probability of the site being a variant is performed, there is no estimate of the selection coefficient (gamma) and its confidence intervals. Please check tutorial for more details!"<<endl;
    cout<<"*************"<<endl;
  }
  int rs=0;
  int RecurSize=Recurrents.size();
  for (rs=0;rs<RecurSize;rs++)
  {
	   int pos=Recurrents[rs].sites;
	   div_codon_consensus[pos]='Q'; //Record recurrent sites as Q
  }
//print CSIMAC main results (gamma, 95% CI gamma) to the output file
  output(N);
  return 1;
}

/***************************************************
Subfunction - open and read the recurrent file
Input: the file name, use the same format of the Recurrent.txt
Output: recurrent list
***************************************************/
int cPRFCluster::GetRecurrentList(string input_f){
	ifstream myfileFn2(input_f.c_str());
	if (!myfileFn2) throw "Error in opening Recurrent File for GetRecurrentList...\n";

	//Read the file; remove the empty lines at the end of the file
	string str;
	while ( myfileFn2.good()) {
		getline(myfileFn2,str);
		// cout<<"Each line: "<<str<<endl;
		int position1 = str.find(" => ");
		if (position1!=std::string::npos)
		{
			//cout<<"The position1: "<<position1<<endl;
		}
		else {cout<<"Error! Failed to find the marker => in the recurrent file for GetRecurrentList!\n";}
		int position2=str.find("	", position1+2);
		//Get the recurrent position and recurrent count
		if (position1<20 and position1>0)
		{
			string e1 = str.substr(0, position1);
			string e2 = str.substr(position1 + 4,
					str.length() - position1 - 4);

			//cout<<e1<<"\t**"<<e2<<"$$\t***"<<endl;
			int pos=CONVERT<int>(e1);
			int count=CONVERT<int>(e2);
			RecurrentSitesCounts tmp_rc(pos, count);
			Recurrents.push_back(tmp_rc);
			//Compact into one site for each codon
		}
	}
	myfileFn2.close();
	/*
   	cout<<"***The size of the recurrent: "<<Recurrents.size()<<endl;
    cout<<"RecurrentSite: "<<Recurrents[0].sites<<"\tCount: ";
    cout<<Recurrents[0].counts<<endl;
	 */
	return 1;
}
/***************************************************
* Function: Count the number for Synonymous or Replacement
* Input Parameter: seq - polymorphism or divergence sequence;start position; end position; symbol - Synonymous (S) or Replacement (R)
* Output: the number of symbols (Synonymous or Replacement)
* Return Value: the number of symbols (Synonymous or Replacement)
***************************************************/
long cPRFCluster::getDifference(string seq, int pos_start, int pos_end, char symbol) {
		long i, n = 0;
		for (i=pos_start; i<=pos_end; i++) {
			if (seq[i]==symbol) n++;
		}
  if (symbol == 'S')
		return n;
	else if (symbol == 'R')
	{
		for (long i = 0; i < Recurrents.size(); i++)
		{
			if (Recurrents[i].sites >= pos_start and Recurrents[i].sites <= pos_end)
			{
				n += (Recurrents[i].counts - 1);
			}
		}
		return n;
	}
	else return -1;
}
/***************************************************
* Function: calculate the likelihood of Binomial Probability i*log(p)+(n-i)*log(1-p); p=i/n; Log likelihood of Bernoulli distribution. BernoulliProb= (i/n)^i*[(n-i)/n]^(n-i); prob=Log(BernoulliProb)=i*log(i/n)+(n-i)*log[(n-i)/n]
* Input Parameter: total number n and occurence i
* Output: probability
* Return Value: probability
***************************************************/
double cPRFCluster::BinomialProb(long n, long i) {
  double prob = 0.0;
  prob += (i==0)?0:i*log(double(i)/n);
  prob += (n-i==0)?0:(n-i)*log(double(n-i)/n);
  return prob;
}
/***************************************************
* Function: Calculate the log of factorial, to prevent the large number.
***************************************************/
double cPRFCluster::factorial(int n) {
	//ZMZ added 04/26/2016
	double ans = 1;
	//ZMZ added 04/25/2016
	if (n<1) { 
	cout<<"\nRecurrent Count in Factorial: "<<n<<endl;
	throw "Error in Factorial with negative values from Recurrent!\n";}
	for (int i = 1; i <=n; i++)
		ans += log(i);
	return ans;
}
/***************************************************
* Function: Calculate loglikelihood for the cluster using the poisson rate
***************************************************/
double cPRFCluster::LogLikelihoodCluster(long cs, long ce, long start, long end, char symbol) {
	//ZMZ added 04/26/2016
	double M = 0;
	long N_cluster_ScaledBack = (ce - cs + 1)*Scale;
	long n = getDifference(div_codon_consensus,cs,ce,symbol);
	double lambda = (double)n/N_cluster_ScaledBack; //poisson rate
	double likelihood;	
	if (n==0) {	//To prevent log(0) case, get the pseudo-lambda
		lambda = (double)1/(N_cluster_ScaledBack+1); //ZMZ 07/18/2016. Fixed bug, added double, or it will be 0. lambda = (double)1/(N_cluster_ScaledBack+1);
		likelihood=lambda * N_cluster_ScaledBack * log(lambda) - N_cluster_ScaledBack * lambda; 
		//cout<<"\nTest LogLikelihoodCluster, lambda:\t"<<lambda<<"\tlog(lambda):\t"<<log(lambda)<<"likelihood:"<<likelihood<<"\t- n:"<<-n<<endl;
		//return likelihood; 
		return 0;
		}		
	if (symbol=='S') {//for synonymous, recurrent is not considered; or the replacement recurrent M will be calculated, and create a bug.
		M=0;
	    likelihood=n * log(lambda) - n - M;
		return likelihood; 
 	}
	for (long i = 0; i < Recurrents.size(); i++)
	{
		if (Recurrents[i].sites >= cs and Recurrents[i].sites <= ce)
		{
			//ZMZ added 04/26/2016
			int RecCount=CONVERT<int>(Recurrents[i].counts);
			//n += (Recurrents[i].counts - 1); calculation moved to getDifference
			M += factorial(RecCount);
			//ZMZ added 04/26/2016
			if (RecCount<2) {
			cout<<"\nLogLikelihoodCluster: Recurrent Count:$"<<Recurrents[i].counts<<"$Converted Recurrent Count:$"<<RecCount<<"$Recurrent Factorial:$"<<M<<endl;
			throw "LogLikelihoodCluster: Error in Recurrent Factorial: not integer Recurrent Count!\n";
			}	
		}
	}
	likelihood=n * log(lambda) - n - M;
	//cout<<"\nTest LogLikelihoodCluster, lambda:\t"<<lambda<<"likelihood:"<<likelihood<<"\tlog(lambda):\t"<<log(lambda)<<"\t- n:"<<-n<<"\t-M:"<<-M<<endl;
	return  likelihood;
}
/***************************************************
* Function: Calculate loglikelihood for the non-clusters using the poisson rate
***************************************************/
double cPRFCluster::LogLikelihoodNonCluster(long cs, long ce, long start, long end, char symbol) {
	//ZMZ added 04/26/2016
	double M = 0;
	long N_noncluster_ScaledBack = ((end-start + 1)-(ce - cs + 1))*Scale;
	long n = getDifference(div_codon_consensus,start,cs-1,symbol) + getDifference(div_codon_consensus,ce+1,end,symbol);
	double lambda = (double)n/N_noncluster_ScaledBack;
	double likelihood;
	if (n==0) {	//To prevent log(0) case, get the pseudo-lambda
		lambda = (double)1/(N_noncluster_ScaledBack+1); 
		likelihood=lambda * N_noncluster_ScaledBack * log(lambda) - N_noncluster_ScaledBack * lambda;
		//return likelihood;
		return 0;
		} 
	if (symbol=='S') {//for synonymous, recurrent is not considered; or the replacement recurrent M will be calculated, and create a bug.
		M=0;
	    likelihood=n * log(lambda) - n - M;
		return likelihood; 
 	}		
	for (long i = 0; i < Recurrents.size(); i++)
	{
		if (Recurrents[i].sites < cs || Recurrents[i].sites > ce)
		{
			//ZMZ added 04/26/2016
			int RecCount=CONVERT<int>(Recurrents[i].counts);
			//n += (Recurrents[i].counts - 1); calculation moved to getDifference
			M += factorial(RecCount);
			//ZMZ added 04/26/2016
			if (RecCount<2) {
			cout<<"\nLogLikelihoodNonCluster: Recurrent Count:$"<<Recurrents[i].counts<<"$Converted Recurrent Count:$"<<RecCount<<"$Recurrent Factorial:$"<<M<<endl;
			throw "nLogLikelihoodNonCluster: Error in Recurrent Factorial: not integer Recurrent Count!\n";
			}	
		}
	}	
	likelihood=n * log(lambda) - n - M;
	//cout<<"\nTest LogLikelihoodNonCluster, lambda:\t"<<lambda<<"likelihood:"<<likelihood<<"\tlog(lambda):\t"<<log(lambda)<<"\t- n:"<<-n<<"\t-M:"<<-M<<endl;
	return likelihood; 
}
/***************************************************
* Function: Determine the hot and cold spots by calculating the percentage of symbol counts in the cluster (nc/cent_len) and non-cluster regions (nw-nc)/non_cent_len.
* Input Parameter: sequence start and end as pos_start, pos_end;cluster start and end as cs and ce; probablity of cluster and non-cluster region as pc and p0; number of symbol (Synonymous or Replacement) counts in the whole sequence and in the cluster only (nw and nc)
* Output: percentage of symbols in the cluster and non-cluster regions
* Return Value:
***************************************************/
double cPRFCluster::getp0pc_MK(int pos_start, int pos_end, int cs, int ce, float &p0, float &pc, int nw, int nc) {
  int non_cent_len=pos_end-pos_start-ce+cs; //the length for non-cluster sequence.
  int non_cent_len_ScaledBack=non_cent_len*Scale;  //Scale the length for non-cluster sequence back   
  if(cs==pos_start && ce==pos_end){
    p0=0.0; //rate of variants in the non-cluster region is 0 if the cluster part is the whole sequence and the non-cluster sequence is empty.
  }else{
    //nw means the total number of variants in the whole sequence; nc means the number of variants in the cluster region.
	  p0=(float)(nw-nc)/non_cent_len_ScaledBack; // rate of variants in the non-cluster region
  }
  //pc means the cluster part; nc means number of symbols (Synonymous or Replacement) in the cluster region
  int cent_len=ce-cs+1; // the length for the cluster
  int cent_len_ScaledBack=cent_len*Scale;  //Scale the length for the cluster back 
  pc=(float)nc/cent_len_ScaledBack; // rate of variants in the cluster region
  return 1;
}
/***************************************************
* Function: Find cluster and calculate likelihood of sites being variant under different models AIC, BIC, AICc
* Input Parameter: start position, end position, symbol - synonymous or replacement,?SiteModels *pointer?
* Output: vectors vec_AllModels; vec_SelectedModels
* Return Value:
***************************************************/
int cPRFCluster::ClusterSubSeq(int pos_start, int pos_end, char symbol, struct SiteModels *pointer) {
	time_t time_start1 = time(NULL); // Record the start time
	long N = pos_end - pos_start + 1; // total region length
	long N_ScaledBack=N*Scale; // Scale the gene length back 
	long symbol_n = 0; // declare the counts of variant sites
	symbol_n=getDifference(div_codon_consensus, pos_start, pos_end, symbol); // Get the counts of variant sites
	if (N==0 || N==1) 
	{
		cout<<"No cluster; quit ClusterSubSeq for the region "<<pos_start<<"\tto\t"<<pos_end<<"\tTotalVariant\t"<<symbol_n<<endl;
		return 1; // return if the sequence length is 0/1 or the count of variant sites is 0	
	}

	if (symbol_n==0 && flag_found_dr==0) //when n=0 and no cluster found ever, quit cluster; else keep all cluster models for sub-regions.
	{
		cout<<"No cluster; quit ClusterSubSeq for the gene "<<pos_start<<"\tto\t"<<pos_end<<"\tTotalVariant\t"<<symbol_n<<endl;
		return 1; // return if the sequence length is 0/1 or the count of variant sites is 0	
	}	
	
	double InL0 = LogLikelihoodCluster(0, N-1, 0, N-1,symbol); // the background likelihood for the whole gene/region
	double InL_cluster_max, InL_noncluster_max; // for cases with 0 variant log(0)
	InL_cluster_max=InL_noncluster_max=InL0;
	double AIC0,AICc0,BIC0,cri0;
	cri0=AIC0=AICc0=BIC0=-2*InL0; // the overall weight for the whole region
	double InL,AIC,AICc,BIC;
	InL = InL0;
	AIC = AICc = BIC = AIC0; // Initialization based on the whole region, parameter is 0 for the whole region with the cluster in the whole region.

	double min_cri = 10000000; // a big value for an initial of the minimal criteria to keep the best model, the smaller cri, the better the model.
	int para=2; // the number of parameters, cs and ce
	long cs, ce, cs_max, ce_max; //cluster start and end positions
	double lambda_0_max,lambda_c_max;
	lambda_0_max=0;
	lambda_c_max=symbol_n/N_ScaledBack;

	int found=0; // means the absence of the cluster; found =1, found the lowest AIC/BIC, as the best model
	vec_AllModels.clear(); // empty the vec_AllModels for the new region
	//ZMZ debugging 04/27/2016
	//cout<<"Model Information\npos_start\tpos_end\tcs\tce\tp0\tpc\tsymbol_n\tsymbol_cn\tInL_tmp\tInL_tmp_cluster\tInL_tmp_noncluster\tcri\tcri0\tLnL0"<<endl;

	cout<<"Within the clustering subfunction ClusterSubSeq, start of the region: "<<pos_start<<"\tEnd: "<<pos_end<<"\tTotalVariant\t"<<symbol_n<<endl;
	//Slide window across the whole sequence for the cluster start and end position to find the cluster with the window size as 1.
	for (cs=pos_start; cs<=pos_end; cs+=1) {
		for(ce=cs; ce<=pos_end; ce+=1) {
			para=2;
			if (cs==pos_start && ce==pos_end) para = 0;//the number of parameters is zero, since both cs and ce are known
			else if(ce==pos_end || cs==pos_start) para=1; //the number of parameters is one, since one of cs and ce is known

			long symbol_cn = 0; //declare the counts of variant sites in the cluster region only.
			symbol_cn = getDifference(div_codon_consensus, cs, ce, symbol);//Get the counts of variant sites in the cluster region only.
			long symbol_ncn = 0; //declare the counts of variant sites in the non-cluster region only.
			symbol_ncn = getDifference(div_codon_consensus, pos_start, cs, symbol)+ getDifference(div_codon_consensus, ce, pos_end, symbol);//Get the counts of variant sites in the non-cluster region only.
			long total_variant=symbol_cn + symbol_ncn;

			//Calculate log likelihood of the sites being variant in the cluster & Non-cluster			
			double InL_tmp_cluster = LogLikelihoodCluster(cs, ce, pos_start, pos_end,symbol); 
			double InL_tmp_noncluster = LogLikelihoodNonCluster(cs, ce, pos_start, pos_end,symbol);
			//if (symbol_ncn==0) {InL_tmp_noncluster=InL0;}	
			//if (symbol_cn==0) {InL_tmp_cluster=InL0;}
				
			double InL_tmp = InL_tmp_cluster + InL_tmp_noncluster;
			double AIC_tmp  = -2*InL_tmp + 2*para; // AIC=-2ln(L)+2k
			double AIC_weight=AIC_tmp; // AIC weight
			double AICc_tmp = AIC_tmp;
			if (N_ScaledBack-para-1>0.0) AICc_tmp += 2*para*(para+1)/(N_ScaledBack-para-1); //calculate AICc=AIC+2k(k+1)/(l-k-1)
			else AICc_tmp = 2*AIC_tmp;
			double BIC_tmp = -2*InL_tmp + para*log(double(N_ScaledBack)); //calculate BIC=-2ln(L)+kln(l)
			double cri; //criteria for the model selection

			//AIC, default
			if (criterion_type==0){
				cri = AIC_tmp;
				cri0=AIC; //update with the best model
			}
			//BIC
			else if (criterion_type==1){
				cri = BIC_tmp;
				cri0=BIC; //update with the best model
			}
			//AICc
			else if (criterion_type==2){
				cri = AICc_tmp;
				cri0=AICc; //update with the best model
			}

			if(cri < min_cri) min_cri=cri;//Get the minimal value of the selected model or criteria, which corresponds to the best model.
			float p0, pc; // p0 and pc are the percentage of symbols in the non-cluster and cluster regions; they are used to determined if the cluster is under a hot (p0<pc) or cold (p0>pc) spot.
			getp0pc_MK(pos_start, pos_end, cs, ce, p0, pc, symbol_n, symbol_cn);
			CandidateModels tmp_CM(AIC_weight, cri, pos_start, pos_end, cs, ce, p0, pc,InL_tmp);// each model contains all parameters
			vec_AllModels.push_back(tmp_CM); // all models
			/*
			if (cs==170 and ce==289){
				cout<<"Cluster region 170 to 289."<<endl;
				cout<<"Position Start: "<<pos_start<<"\tEnd: "<<pos_end<<"\tCluster start: "<<cs<<"\tend: "<<ce<<"\tp0: "<<p0<<"\tpc: "<<pc<<"\tvariant#Total\t"<<symbol_n<<"\tVariant#InCluster:\t"<<symbol_cn<<"\tInL0:\t"<<InL0<<"\tInL_tmp:\t"<<InL_tmp<<"\tInL_tmp_cluster:\t"<<InL_tmp_cluster<<"\tInL_tmp_noncluster:\t"<<InL_tmp_noncluster<<endl;		
			}*/

			//Evaluate the cluster by the criterion, found the best cluster model, the first cri0 is Null model cs=pos_start, ce=pos_end.
			if (cri < cri0 && ce-cs > 1 && symbol_n > 1) {//add the condition to exclude the cluster of one site, and no cluster for only one variant
				if ((cs-pos_start>1 && pos_end-ce>1) ||
						(cs==pos_start && pos_end-ce>1) ||
						(cs-pos_start>1 && pos_end==ce))
				{
					found = 1; // found=1 means the presence of the cluster
					//cout<<"\n***Found a cluster (cri < cri0), cri: "<<cri<<"\tcri0: "<<cri0<<endl;
					//ZMZ debugging 07/14/2016
					//cout<<"Position Start: "<<pos_start<<"\tEnd: "<<pos_end<<"\tCluster start: "<<cs<<"\tCluster end: "<<ce<<"\tp0: "<<p0<<"\tpc: "<<pc<<"\tvariant#Total\t"<<symbol_n<<"\tVariant#InCluster:\t"<<symbol_cn<<"\tInL0:\t"<<InL0<<"\tInL_tmp:\t"<<InL_tmp<<"\tInL_tmp_cluster:\t"<<InL_tmp_cluster<<"\tInL_tmp_noncluster:\t"<<InL_tmp_noncluster<<endl;
					cs_max = cs;
					ce_max = ce;
					lambda_0_max=p0; 
					lambda_c_max=pc;
					//Updates the cri0 with the current cluster model; overall find the best model in the region pos_start to pos_end
					InL = InL_tmp;
					AIC = AIC_tmp;
					AICc = AICc_tmp;
					BIC = BIC_tmp;
				}
			}//end of if (cri < cri0)
		} // end of the inside for loop, sliding windows of the end of the cluster within the region
	}//end of the outside for loop, sliding windows of the start of the cluster within the region
	cout<<"For the region from "<<pos_start<<" to "<<pos_end<<", the total number of models in the subfunction ClusterSubSeq: "<<vec_AllModels.size()<<endl;
	// if no cluster is found, vec_SelectedModels keeps the null model cs=pos_start, ce=pos_end, and then do EachSiteModels.
	if (found==0){
		//If no cluster is found in the whole gene, then keep the null model.
		if((symbol=='S' && flag_found_ds==0) || (symbol=='R' && flag_found_dr==0)){
			double p_tmp=(double)symbol_n/(double)N_ScaledBack;
			CandidateModels nullmodel(0, N-1, 0, N-1, p_tmp,p_tmp, InL0, InL, AIC0, AIC, AICc0, AICc, BIC0, BIC);
			vec_SelectedModels.push_back(nullmodel);
			cout<<"No cluster is found for the entire gene; Model averaging of models for each site. exit current ClusterSubSeq."<<endl<<endl;
			EachSiteModels(0, N-1, 0, N-1, p_tmp,p_tmp, min_cri,pointer); //Keep the null model for each site in the vector pointer if no cluster found in the whole gene
			return 1;
		}
		else {
			EachSiteModels(pos_start, pos_end, cs_max, ce_max, lambda_0_max, lambda_c_max, min_cri,pointer); //keep all models in the sub-region even no cluster is found
			cout << "No cluster is found for the subregion: exit current ClusterSubSeq."<<endl;
			return 1;
		}
	}
	// if the cluster is found, keep the cluster in vec_SelectedModels, and then do EachSiteModels.
	else{
		if(symbol=='S') flag_found_ds++;//divergence synonymous
		if(symbol=='R') flag_found_dr++;//divergence replacement
	}

	//in the case that the cluster is present, Select models (vec_SelectedModels) keep the best model cs, ce; Model average.
	CandidateModels selectedmodel(pos_start, pos_end, cs_max, ce_max, lambda_0_max, lambda_c_max, InL0, InL, AIC0, AIC, AICc0, AICc, BIC0, BIC);
	vec_SelectedModels.push_back(selectedmodel); //Accumulate the best clustering model or null model for each region, only one selected for each region in each ClusterSubSeq.
	cout<<"Found a cluster. The selected model: "<<pos_start<<"-"<<pos_end<<"\tCluster: "<<cs_max<<"-"<<ce_max<<endl;
	cout<<"Total number of total selected models in ClusterSubSeq: "<<vec_SelectedModels.size()<<endl;
	cout << "Model averaging of models for each site, to get the probability of the site being variants."<<endl;
	EachSiteModels(pos_start, pos_end, cs_max, ce_max, lambda_0_max, lambda_c_max, min_cri,pointer); //Keep all models for each site in the vector pointer
	

	/* Divide and Conquer: do ClusterSubSeq for the three sub-sequences (pos_start to cs, ce to pos_end, and cs to ce) for the best models cs, ce*/
	if (ce_max!=pos_end || cs_max!=pos_start) {
		if (cs_max>pos_start+1) {
			cout<<"Divide and Conquer: do ClusterSubSeq for regions under cs_max>pos_start+1: "<<pos_start<<" — "<<(cs_max-1)<<endl;
			ClusterSubSeq(pos_start, cs_max-1,symbol,pointer);}
		if (ce_max<pos_end-1) {
			cout<<"Divide and Conquer: do ClusterSubSeq for regions under ce_max<pos_end-1: "<<ce_max+1<<" — "<<pos_end<<endl;
			ClusterSubSeq(ce_max+1, pos_end,symbol,pointer);}
		if (cs_max<ce_max-1) {
			cout<<"Divide and Conquer: do ClusterSubSeq for regions under cs_max<ce_max-1: "<<cs_max<<" — "<<ce_max<<endl;
			ClusterSubSeq(cs_max, ce_max, symbol,pointer);}
	}
	cout<<"Finish ClusterSubSeq for the region from the start position "<<pos_start<<" to the end position "<<pos_end<<endl;
	return 1;
}
/***************************************************
 * Function: Model Average for each site in the region pos_start and pos_end
 * Input Parameter: region pos_start and pos_end, the best model (cs_max, ce_max), the likelihood for the best model (lambda_0_max, lambda_c_max), best model AIC (min_cri)
 * Output: model averaged rate for each site in the region; SiteModels *pointer vector (to calculate CI_MA) contains the site, all models for the site, each model with AIC weight, Site rate, LogLikelihood
 * Return Value: 1
 ***************************************************/
int cPRFCluster::EachSiteModels(long pos_start, long pos_end, long cs_max, long ce_max, double lambda_0_max, double lambda_c_max, double min_cri, struct SiteModels *pointer){
	long i, j;
	
	/* Under model selection, assign lambda, the poisson rate for each position based on the best cluster model in the region pos_start and pos_end. */
	for (i=pos_start; i<=pos_end; i++) {
		vec_MS_rate[i] = lambda_0_max; //For site not in the cluster, the MS rate is the non-cluster rate
	}
	for (i=cs_max; i<=ce_max; i++) {
		vec_MS_rate[i] = lambda_c_max; //For sites in the cluster (between cs and ce), the MS rate is the cluster rate.
	}

	/* Model Averaging rate based on all models for the region pos_start and pos_end*/
	double all_weight = 0.0;
	//Get the AIC weight for each model i
	for (i=0; i<vec_AllModels.size(); i++) { //Calculate the summed weight all_weight
		vec_AllModels[i].AICweight = exp(-0.5*(vec_AllModels[i].AICweight - min_cri));
		all_weight += vec_AllModels[i].AICweight;
	}
	for (i=0; i<vec_AllModels.size(); i++) { //Calculate the weight for each model i
		vec_AllModels[i].AICweight = vec_AllModels[i].AICweight/all_weight;
	}

	//for each site i, get the model-averaged poisson rate and and all models j in pointer
	for (i=pos_start; i<=pos_end; i++) {
		double rate = 0.0;
		vector<CI> CIs;
		CIs.clear();

		//Model average across all models, and get site i and model j specific parameters (AICweight, poisson rate, loglikelihood). 
		for (j=0; j<vec_AllModels.size(); j++) {
			double site_i_model_j_AICweight=0.0;
			site_i_model_j_AICweight = vec_AllModels[j].AICweight; // AIC weight
			double site_i_model_j_rate=0.0;
			double site_i_model_j_Loglikelihood=vec_AllModels[j].LogLikelihood;
			if (i<vec_AllModels[j].cs || i>vec_AllModels[j].ce) {
				site_i_model_j_rate = vec_AllModels[j].p0; //site rate for the specific model if the site is not in the cluster in that model
			}
			else {
				site_i_model_j_rate = vec_AllModels[j].pc; //site rate for the specific model if the site is in the cluster in that model
			}
			//cout<<"Model\t"<<j<<"\tAIC weight:\t"<<site_AICweight<<"\tLoglikelihood\t"<<site_likelihood<<"\tSiteRate\t"<<site_rate<<endl;
			CI tmp_cis(site_i_model_j_AICweight, site_i_model_j_rate,site_i_model_j_Loglikelihood); //for each model, keep the AIC weight, rate, logLikelihood.
			CIs.push_back(tmp_cis); // keep all models for the site
			rate += site_i_model_j_rate*site_i_model_j_AICweight; //model averaged rate across all models
		}//end of the inside for loop, for each model under one position.
		vec_MA_rate[i] = rate; // assign the model averaged rate for the site
		// pointer vector contains the site, all models for the site, each model with AIC weight, site rate, LogLikelihood
		//pointer[i].sms.clear();// silent pointer[i].sms.clear() in EachSiteModels, since cumulative models are used in ClusterSubSeq for each site.
		pointer[i].sms=CIs;//Keep weight, poisson rate, and likelihood for each model sms for each site i
		pointer[i].pos=i; //Keep the position of each site.
		CIs.clear();
	}//end of the outside for loop, for each position
	return 1;
}
/***************************************************
* Function: Get Confidence Interval for Model Average
* Input Parameter:
* Output: vectors vec_lower_rate, vec_upper_rate
* Return Value: none
***************************************************/
int cPRFCluster::CI_MA(struct SiteModels *pointer,long N){
  for(long i=0;i<N;i++){
    //probability and weight for all possible models
    vector<CI> CIs;
    CIs.clear();
    CIs=pointer[i].sms;
    double lower=0.0;
    double upper=0.0;
    while (lower<quantile_for_CI || upper<quantile_for_CI) {
      long min_pos=0, max_pos=0;
      long j = 0;
      //flag_lower==0, max_pos is in front of min_pos
      //flag_lower==1, max_pos is behind of min_pos and min_pos will be deleted.
      int flag_lower=0;
      while (j<CIs.size()) {
        if(CIs[j].p<CIs[min_pos].p) min_pos = j;
        if(CIs[j].p>CIs[max_pos].p) max_pos = j;
        j++;
      }

      if(lower<quantile_for_CI){
        if(max_pos>min_pos){
          flag_lower=1;
        }

        lower += CIs[min_pos].weight; // add up weight till lower=0.025, by removing all lower CI
        //Get the lower CI probability rate at 0.025
        vec_lower_rate[i] = CIs[min_pos].p;
        CIs.erase(CIs.begin() + min_pos); //truncate the list
      }

      if(upper<quantile_for_CI){
        if(flag_lower==1 && max_pos!=0){
          max_pos--;
        }
        upper += CIs[max_pos].weight; // add up weight till upper=0.025
        //Get the upper CI probability rate at 0.975
        vec_upper_rate[i] = CIs[max_pos].p;
        CIs.erase(CIs.begin() + max_pos);
      }
      flag_lower=0;

      //For some sites,one weight is significant hight than others, so lower=upper
      if(CIs.size()==0){
        if(upper > quantile_for_CI){
          vec_lower_rate[i]=vec_upper_rate[i];
        }else if(lower > quantile_for_CI){
          vec_upper_rate[i]=vec_lower_rate[i];
        }
        break;
      }
      if(CIs.size()==1){
        if(lower<quantile_for_CI && upper<quantile_for_CI){
          vec_lower_rate[i]=CIs[0].p;
          vec_upper_rate[i]=CIs[0].p;

          break;
        }
      }    
    }             
    CIs.clear();
  }
  return 1;
}

/***************************************************
* Function: Calculate replacement mutation rate based on synonymous rate and the estimated ratio of Replacement and Synonymous (N/S) based on human reference sequence for a gene
* Input Parameter: N/S ratio; synonymous mutation rate
* Output: replacement mutation rate
* Return Value: replacement mutation rate
***************************************************/
double cPRFCluster::ReplacementRate(double ratio_NS, double syn_rate){
  //double ur=ratio_NS*syn_rate;//calculate replacement rate. This is wrong
  double ur=syn_rate/ratio_NS;//calculate replacement rate
  //cout<<"In ReplacementRate, Nonsysnonymous and Synonymous ratio ratio_NS: "<<ratio_NS<<endl<<"Synonymous rate syn_rate: "<<syn_rate<<endl<<"Replacement mutation rate ratio_NS*syn_rate: "<<ur<<endl;
  return(ur);
}

/***************************************************
* CancerSynonymousRate: Use the silent rate from MutSigCV, corrected silent mutation rate with gene expression,replicating time,and chromosomal position.
* Function: Calculate cancer Divergence Silent mutation rate based on Formula us(c)= NumberOfSilentMutation/(tumor_num*GeneLength)
* Input Parameter: the sample number of cancer divergence, and the sequence length
* Output: cancer divergence silent mutation rate
* Return Value: cancer divergence silent mutation rate
***************************************************/
double cPRFCluster::CancerSynonymousRate(long tumor_num, long N){
  double ds=0.0; //count of synonymous sites in the cancer divergence sequence
  ds=getDifference(div_codon_consensus,0,N-1,'S'); // Count cancer divergence sysnonymous sites

  double us_c=0.0;
  int N_ScaledBack=N*Scale; //Scale the gene length back 

  //Use the input MutSigCV mutation rate as the priority, if no input, calculate silent rate using silent mutations in the gene
  if (SilentRate_flag==1)
  {
	  us_c=SilentRate;
  }
  else
  {
	  if (ds!=0.0)
	  {
		  // Calculate cancer divergence silent mutation rate based on Formula us(c)= NumberOfSilentMutation/(tumor_num*GeneLength)
		  cout<<"Warning: No input silent mutation rate.Calculate cancer divergence silent mutation rate based on Formula us(c)= NumberOfSilentMutation/(tumor_num*GeneLength)\n";
		  us_c=ds/(N_ScaledBack*tumor_num);
	  }
	  else
	  {
		  cout<<endl<<endl<<"Error. The number of synonymous mutation is 0. The genomic average is required (-ga value) as an input."<<endl<<"See [./CSI-MAC -h] for more information."<<endl<<endl<<endl;
		  throw 1;
	  }
  }
  cout<<"Cancer Divergence synonymous mutation rate: "<<us_c<<endl;
  return(us_c);
}

/***************************************************
* Function: Exact algorithm estimating gamma and confidence intervals for parameter r at each site
* Input Parameter: dr - all models for all sites, each model contains (site_rate, site_AICweight, LogLikelihood); Gene length - N
* Output:
* Return Value: vec_r_c, vec_upper_r_c,vec_lower_r_c
***************************************************/
int cPRFCluster::CIr_exact(struct SiteModels *dr, long N){
  cout<<"Executing the exact algorithm for calculation of the model-averaged estimate of gamma and CIs for gamma"<<endl;
  vector<thread> exact_threads;
  int index=0;
  ostringstream* thread_outputs = new ostringstream[NUM_CPU];

  while(index<N)
  {
    int num_threads=0;
    for (int i=0;(i<NUM_CPU && index<N);i++)
    {
     exact_threads.push_back(thread(&cPRFCluster::CIr_exact_threaded, 
      this, dr, N, index,&thread_outputs[i]));
     index++;
     num_threads++;
    }
   for (long i=0;i<num_threads;i++){
    exact_threads[i].join();
    cout << thread_outputs[i].str();
    thread_outputs[i].str("");
    thread_outputs[i].clear();
    }
    exact_threads.clear();
  }
  return 1;
}

void cPRFCluster::CIr_exact_threaded(struct SiteModels *dr, long N, long i, ostringstream* myout){
	long jj;
	double min_weight_c=1000000;
	double new_r_c=0;
	//Get all the models
	vector<rModels> vec_rModels_c_indiv;
	for(long k=0;k<dr[i].sms.size();k++){
		if(dr[i].sms[k].p==0.0){
			//continue;
			new_r_c=0; //ZMZ 07/15/2016 debug 
		}
		else {
			//dr[i].sms[k].p is the poisson rate of divergence for the site i in the specific model k
			// JT: I think that .p should be .lambda as it is now a rate not a probability
			new_r_c=CIs_rc_PRF(ucr,dr[i].sms[k].p,tumor_num); // calculate gamma for this site under this model		
		}
		double weight_tmp_c=dr[i].sms[k].weight;// Assigns the weight for this gamma at site i under model k as the model's weight.
		if(weight_tmp_c<min_weight_c) min_weight_c=weight_tmp_c;
		rModels tmp_rm_c(weight_tmp_c,new_r_c);
		vec_rModels_c_indiv.push_back(tmp_rm_c);
	}
	//Calculate the CIs for r using exact algorithm
	if(vec_rModels_c_indiv.size()==0){
		vec_lower_r_c[i]=0;
		vec_upper_r_c[i]=0;
		vec_r_c[i]=0;
	}
	else{
		CI_UpLow_rc(i,min_weight_c,vec_rModels_c_indiv,myout); // Get the upper and lower boundaries of gamma for site i by exact algorithm; and model averaged gamma
	}
	vec_rModels_c_indiv.clear();	
}

/***************************************************
* Function: Stochastic algorithm (or exact for models less than 100) estimating confidence intervals for parameter gamma (r) at each site
* Input Parameter: dr - all models for all sites, each model contains site_mutation_rate, site_AICweight, and LogLikelihood; Gene length - N
* Output: 10000 models with gamma and gamma weight; model averaged gamma and 95% CI gamma
* Return Value: vec_r_c, vec_upper_r_c,vec_lower_r_c
***************************************************/
int cPRFCluster::CIr_stochastic( struct SiteModels *dr, long N){
  time_t time_start1 = time(NULL); // Record the start time
  cout<<"******CIr_stochastic******"<<endl;
  vector<thread> stochastic_threads;
  int index=0;
  ostringstream* thread_outputs = new ostringstream[NUM_CPU];

  while(index<N)
  {
    int num_threads=0;
    for (int i=0;(i<NUM_CPU && index<N);i++)
    {
     stochastic_threads.push_back(thread(&cPRFCluster::CIr_stochastic_threaded, 
      this, dr, N, index, time_start1, &thread_outputs[i]));
     index++;
     num_threads++;
    }
   for (long i=0;i<num_threads;i++){
    stochastic_threads[i].join();
    cout << thread_outputs[i].str();
    thread_outputs[i].str("");
    thread_outputs[i].clear();
    }
    stochastic_threads.clear();
  }
  return 1;
}

void cPRFCluster::CIr_stochastic_threaded(struct SiteModels *dr, long N, long i, time_t time_start1, ostringstream* myout){
	//cout<<" Sequence length: "<< N<<" Estimated time: "<<12*N/1000<<" minutes (Time rate: 12 minutes per 1000 sites per 100k models)"<<endl;
	//Estimate gamma CI stochastically for each site
	int parameters=2;
	long jj;
	int RecurrentCount=1;
	//Go over each site of the gene, to calculate gamma and 95% CI gamma
		
    vector<rModels> vec_rModels_c_indiv; // clear the gamma model for the previous site
		//Each time get one dr, estimate r
		long end_num=N_Random; //N_Random=100000
		double min_weight_c=1000000;
		double new_r_c=0; 
		//screen output for every 1000 sites, to monitor the speed of the running
		if (i%1000==0) {
			*myout<<"Start CIr_stochastic at Position "<<i+1;
			*myout<<" (Time elapsed: ";
			time_t t2 = time(NULL)-time_start1;
			int h=t2/3600, m=(t2%3600)/60, s=t2-(t2/60)*60;
			if(h)  *myout<<h<<":"<<m<<":"<<s<<")"<<endl;
			else   *myout<<m<<":"<<s<<")"<<endl;
		}

		//If the number of models for the site is no more than the cutoff N_Random, Use exact CI for gamma.
		if(dr[i].sms.size()<=N_Random){
			*myout<<"Site: "<<i<<"\tModel size: "<<dr[i].sms.size();
			*myout<<"\tExact algorithm is used since the model size is smaller than defined."<<endl;
			////CIr_exact(dr, N);
			//Get all the models
			for(long k=0;k<dr[i].sms.size();k++){
				if(dr[i].sms[k].p==0.0){
					//continue;
					new_r_c=0;
				}
				else {
					// dr[i].sms[k].p is the poisson rate for the site i, the specific model k
					new_r_c=CIs_rc_PRF(ucr,dr[i].sms[k].p,tumor_num); // calculate gamma based on likelihood				
				}
				double weight_tmp_c=dr[i].sms[k].weight;//Same as the MACML AIC weight, don't need to recalculate, just need to re-gather all selected models, and re-calculate the weight of selected models
				if(weight_tmp_c<min_weight_c) min_weight_c=weight_tmp_c;
				rModels tmp_rm_c(weight_tmp_c,new_r_c);
				vec_rModels_c_indiv.push_back(tmp_rm_c);
			}

			//Calculate the CIs for r using exact algorithm
			if(vec_rModels_c_indiv.size()==0){//empty, no gamma values, output -199, NULL
				vec_lower_r_c[i]=-199;
				vec_upper_r_c[i]=-199;
				vec_r_c[i]=-199;
			}
			else{
				CI_UpLow_rc(i,min_weight_c,vec_rModels_c_indiv,myout); // Get the upper and lower boundaries of gamma for site i by exact algorithm; and model averaged gamma
			}
			vec_rModels_c_indiv.clear();
		}
		// if the number of models at one site is larger than the cutoff N_Random, then sort assigned number of random models by AICweight, and compute gamma and stochastic CI for gamma
		else{
			//Sort all different models by weight and probability for a particular site (dr - models for divergence replacement at site i)
			//Rank weight descendingly, if same weight, rank site_rate decreasingly.
			sort(dr[i].sms.begin(), dr[i].sms.end(), more_than_CI());
			////cout<<"Weight for the first model:\t"<<dr[i].sms[0].weight<<endl;
			end_num=N_Random; //N_Random=100000
			//If all the selected models from dr are zero, flag_mutation=0, won't calculate CIs for r.
			int flag_mutation=1;
			int repeat_num=0;
			const vector<double> &pws = RandomModel_NumFastInit(dr,i);
			int TotalSampleNum=N_Random*1000;
			int TotalModel = dr[i].sms.size();
			if (TotalSampleNum>TotalModel) {TotalSampleNum=TotalModel;}//choose the smaller one of N_Random*1000 and Total Model to go through all models, as the cutoff model number when dr.p=0.
			int jj=0; //count all the models sampled, including both p=0 and p!=0
			double new_r_c=-199;
			double weight_tmp_c;
			for(long j=1;j<=end_num;j++){
				//Get a random model by choosing a value from 0 to 1, and select the clustering model for p with the weight above the random value
				//long model_dr=RandomModel_Num(dr,i);
				long model_dr=RandomModel_NumFast(pws); // create a partial summed weight models first, then find the model with the weight according to the random number by NumberOfModels/2
				//discard the model with 0 prob, and don't calculate gamma for it
				jj++; //count all the models sampled, including both p=0 and p!=0
				if(dr[i].sms[model_dr].p==0.0){
					//j--; //repeat getting an alternative model_dr if p=0
					/*
					if(repeat_num==N_Random){
						cout<<"For Site "<<i<<", all assigned selected models indicate no chance to mutate!"; // Every randomly extracted model (all N_Random) gives p=0
						cout<<"\tThe Number of Models with dr.p=0: "<<repeat_num<<";\tThe total Number of Models sampled: "<<TotalSampleNum<<endl;
						repeat_num++;
					}
					*/
					if (jj == TotalSampleNum) //when N_Random samples all have dr.p=0, continue to find random samples till TotalSampleNum
					{
						if (repeat_num == TotalSampleNum) {
							flag_mutation = 0; // only when all models sampled have p=0, the site has no chance to mutation, giving the N-INF for gamma
							*myout << "For Site " << i
									<< ", all sampled models indicate no chance to mutate!"; // Every randomly extracted model (all N_Random) gives p=0
							*myout << "\tThe Number of Models with dr.p=0: "
									<< repeat_num
									<< ";\tThe total Number of Models sampled: "
									<< TotalSampleNum << endl;
						} else {
							flag_mutation = 1; //if there are models p!=0, calculate gamma.
						}
						break;
					} else {
						repeat_num++; //record the number of getting Random model with p=0
					}
					//continue;
					new_r_c=0;
				}
				//keep all non-zero probability models by model specific site rate p and AIC weight, with the assigned number of models as N_Random;
				//calculate gamma and gamma weight based on the p and its weight for all N_Random number of models
				else{
					////RecurrentCount=1; //Not consider the recurrent site at this point
					double site_rate=dr[i].sms[model_dr].p;
					new_r_c=CIs_rc_PRF(ucr,site_rate,tumor_num); // calculate gamma based on likelihood
					}					
				weight_tmp_c=dr[i].sms[model_dr].weight;//Same as the MACML AIC weight, don't need to recalculate, just need to re-gather all selected models, and re-calculate the weight of selected models
				if(weight_tmp_c<min_weight_c) min_weight_c=weight_tmp_c;
				rModels tmp_rm_c(weight_tmp_c,new_r_c);
				vec_rModels_c_indiv.push_back(tmp_rm_c); //keep all calculated gamma and gamma weight for all N_Random number of models
				//cout<<"Model ID:\t"<<j<<"\tSiteRate:\t"<<recur_site_rate<<"\tGamma:\t"<<new_r_c<<"\tAICweight:\t"<<weight_tmp_c<<endl;
				////cout<<"Model ID:\t"<<j<<"\tGamma:\t"<<new_r_c<<"\tAICweight:\t"<<weight_tmp_c<<"\tdr Loglikelihood\t"<<dr[i].sms[model_dr].LogLikelihood<<endl;
				//Debug:check if the models selected are random
				/*
				 if (i==0)
					{cout<<"Model ID:\t"<<j<<"\tGamma:\t"<<new_r_c<<"\tAICweight:\t"<<weight_tmp_c<<"\tdr Loglikelihood\t"<<dr[i].sms[model_dr].LogLikelihood<<"\tRate:\t"<<dr[i].sms[model_dr].p<<endl;
					}
					*/
				
			}
			//Find CI for r using stochastic algorithm
			if(flag_mutation==1){
				// Get gamma, and Upper and Lower CI for gamma using stochastic algorithm
				CI_UpLow_rc(i,min_weight_c,vec_rModels_c_indiv,myout);
			}
			else{
				vec_lower_r_c[i]=-299;
				vec_upper_r_c[i]=-299;
				vec_r_c[i]=-299;
			}
			vec_rModels_c.clear();
		}	
}

/***************************************************
* Function: Get a random model by choosing a value from 0 to 1, and select the clustering model for p with the weight above the random value;
* Function: create a partial summed weight models first, then find the model with the weight according to the random number by NumberOfModels/2
* Input Parameter: SiteModels; site
* Output:model_num
* Return Value: model_num
***************************************************/
vector<double> cPRFCluster::RandomModel_NumFastInit(struct SiteModels *p,long site) {
  int pss = p[site].sms.size();
  CI *cip = &(p[site].sms[0]);
  vector<double> pWeightSums;
  pWeightSums.push_back(cip->weight);
  ++cip;
  for(long i = 1; i < pss; ++i, ++cip) {
    if (cip->weight == 0.0) cerr << "null weight " << site << " " << i << endl;
    pWeightSums.push_back(pWeightSums[i-1] + cip->weight);
  }
  return pWeightSums;
}

long cPRFCluster::RandomModel_NumFast(const vector<double> &pWeightSums){
  double rand_tmp=rand();
  //RAND_MAX is from system (2147483647)
  double rand_num=rand_tmp/RAND_MAX;
  if(rand_num<=pWeightSums[0]){
	  ////cout<<"Warning: the random number is no more than the weight of the first model!"<<endl;
	  ////cout << "RandomNum:\t" << rand_num <<"\tFirstModelWeight:\t"<< pWeightSums[0]<<endl;
      return 0;
  }
  long lo = 0, hi = pWeightSums.size();
  long mid = (lo + hi)/2, nmid;
  while (1) {
    if (pWeightSums[mid] > rand_num)
      hi = mid;
    else
      lo = mid;
    nmid = (lo + hi)/2;
    if (mid == nmid) break;
    mid = nmid;
  }
  //cerr << "Random number: " << rand_num <<" Low: "<< lo << " High:" << hi<< " lowWeightSum: " << pWeightSums[lo] << " highWeightSum: " << pWeightSums[hi] << endl;
  // not exactly the same condition as the original (which i believe is a little broken).
  if(!(hi == 0 || (pWeightSums[hi-1] <= rand_num && rand_num <= pWeightSums[hi]))) {
    cerr << "blew fast model num " << lo << "  " << hi <<  " " << rand_num << ": " << pWeightSums[lo] << " " << pWeightSums[hi-1] << " " << pWeightSums[hi] << endl;
    exit(1);
  }
  ////cout << "RandomNum:\t" << rand_num <<"\tLow:\t"<< lo << "\tHigh:\t" << hi<< "\tlowWeightSum:\t" << pWeightSums[lo] << "\thighWeightSum:\t" << pWeightSums[hi] << endl;
  return(hi);
}

/***************************************************
* Calculate gamma based on the given probability.
* Function: Newton's method to Calculate gamma for cancer divergence by iterations and approximation.
* Input Parameter: cancer replacement rate, ucr; cancer replacement probability calculated by MACML for each site, p_dr; sample number for polymorphism or divergence time, tumor_num
* Output: gamma for cancer
* Return Value: gamma for cancer
***************************************************/
double cPRFCluster::CIs_rc_PRF (double ucr, double dr, long tumor_num){
	//cout<<"******Start CIs_rc_PRF******"<<endl;
	  double j=1.0;
	  double new_r;
	  double min_dx=IR_H;
	  double optimal_r=IR_H;
	  //Estimate r for each site one by one
	  if(dr==0){
		  //Under infinite negative selection
		  //new_r=-299;
		  new_r=0;
		  return (new_r);
	  }else if(ucr==0){
		  new_r=299;
		  return (new_r);
	  }

	  //Calculate the parameter r
	  double tmp=ucr*tumor_num/dr;
	  double df,f;
	  double rtn,dx;

	  bool flag_root=false;
	  for(rtn=IR_L;rtn<=IR_H;rtn+=0.5){
		  if(flag_root==true){
			  break;
		  }
		  for(j=1; j<=JMAX; j++){
			  double tmp0=-2*rtn; //rtn is gamma
			  f=(1-exp(tmp0))/(2*rtn)-tmp;//f(r0) after trapezoidal rule numeric integration of x
			  df=(2*rtn*exp(tmp0)-(1-exp(tmp0)))/(2*rtn*rtn);//f(r0)'after trapezoidal rule numeric integration of x
			  dx=f/df; //f(r0)/f(r0)'
			  rtn -= dx; // iteration for gamma in Newton's method, r1=r0-f(r0)/f(r0)'
			  if(fabs(dx)<fabs(min_dx)){
				  min_dx=dx;
				  optimal_r=rtn;
			  }
			  //criteria to quit the iterations for gamma calculation.
			  if(fabs(dx)<ER){
				  new_r=rtn;
				  flag_root=true;
				  break;
			  }
		  }
	  }
	  if(flag_root==false and min_dx>MinDx){
		  cout<<"CIs_rc_PRF: Gamma NULL. min_dx:\t"<<min_dx<<"\tOptimalGamma:\t"<<optimal_r<<endl;
		  new_r=-199;
	  }
	  if(flag_root==false and min_dx<=MinDx){
		  new_r=optimal_r;
	  }
	  //cout<<"End CIs_rc_PRF: gamma "<<new_r<<endl;
	  return (new_r);
}
/***************************************************
* Function: rank by gamma, create a partial summed weight models first, then find the model with the weight according to the CI
* Function: Use the same one for both Stochastic and Exact to find the gamma CI
* Input Parameter: SiteModels; site; min_weight_c
* Output: lower and upper gamma
* Return Value: int 1
***************************************************/
// this is an expedient lie...
inline bool operator<(const rModels& a, const rModels& b)
{
	return a.r < b.r;
}

int cPRFCluster::CI_UpLow_rc(long site,double min_weight_c, vector<rModels> vec_rModels_c_indiv, ostringstream* myout) {
	double all_weight=0.0;
	long lower_model=-1;
	long upper_model=-1;
	vector<double> rWeightSums_indiv;
    int flag_lower=0;
	double lower=0.0, upper=0.0;

	//Recurrent Site related: Get recurrent weight, recurrentSite gamma, 95% Do_ci_r based on recurrentCount.
	int SiteRecurrentCount=1;
	int jj=0;
	double recur_weight=0.0;
	double recur_r=0.0;
	double recur_upper_r=0.0;
	double recur_lower_r=0.0;

	//ZMZ 04/27/2016 use only poisson rate regional gamma, not recurrent gamma	
	//Find the position in the recurrent list, if find the recurrent site, keep the RecurrentCount
	for (jj=0;jj<Recurrents.size();jj++){
		//cout<<"RecurrentSite:\t"<<Recurrents[jj].sites<<"\tCount: "<<Recurrents[jj].counts<<endl;
		if (Recurrents[jj].sites==site){
			SiteRecurrentCount=Recurrents[jj].counts;
			//cout<<"Site "<<i<<" is in the Recurrent list with RecurrentCount "<<RecurrentCount<<endl;
			//cout<<"SiteRate\t"<<"RecurrentSiteRate\n";
			break;
		}
	}

	//If the site is recurrent site, get the weight for the recurrent site
	if (SiteRecurrentCount>1){
		recur_weight=float(SiteRecurrentCount)/(TotalReplacementSite+TotalRecurrentCount-TotalRecurrentSite); // the weight for the recurrent site
		*myout<<"For Site: "<<site<<"\tRecurrent weight: "<<recur_weight<<"\tSiteRecurrentCount: "<<SiteRecurrentCount<<"\tTotalReplacementSite: "<<TotalReplacementSite<<"\tTotalRecurrentCount: "<<TotalRecurrentCount<<"\tTotalRecurrentSite: "<<TotalRecurrentSite<<endl;
	    //Calculate the gamma and 95% CI gamma for recurrent site, using formula 2r/(1-e^(-2r))=RecurrentNumber/(ReplacementRate*TumorNumber)
	  recurrent_SiteGamma(tumor_num, ucr, SiteRecurrentCount, site, myout); //calculate recurrent site gamma
		double lowCI, upCI;
		LambdaCIs[SiteRecurrentCount-2].count=SiteRecurrentCount;
		lowCI=LambdaCIs[SiteRecurrentCount-2].lowerCIlambda;
		upCI=LambdaCIs[SiteRecurrentCount-2].upperCIlambda;

		recurrent_SiteGammaCI(tumor_num, ucr, lowCI, site,0, myout); // 0 - lowerCI
		recurrent_SiteGammaCI(tumor_num, ucr, upCI, site,1, myout); //1 - upperCI

		recur_r=vec_r_c_r[site];
		recur_upper_r=vec_upper_r_c_r[site];
		recur_lower_r=vec_lower_r_c_r[site];
		*myout<<"Recurrent gamma before correction by weight: "<<recur_r<<"\tlower_r: "<<recur_lower_r<<"\tupper_r: "<<recur_upper_r<<endl;
	}

	sort(vec_rModels_c_indiv.begin(), vec_rModels_c_indiv.end()); // sort by gamma values incrementally
    //debug
	/*
	cout<<"Sorted gamma models:\n"<<"ModelID\tWeight\tGamma\n";
	for (long ii=0; ii<10000;ii++){
		if (ii%100==0){
			cout<<ii<<"\t"<<vec_rModels_c[ii].weight<<"\t"<<vec_rModels_c[ii].r<<endl;
		}
	}
	*/

	//Need to re-calculate the weight, since different set of models are selected, but the same AIC weight can be used as MACML clusters.
	//Calculate weight for each model for the total random picked models and get a summed weight for gamma
	for(long i=0;i<vec_rModels_c_indiv.size();i++){
		vec_rModels_c_indiv[i].weight=exp(-0.5*(vec_rModels_c_indiv[i].weight-min_weight_c));
		all_weight+=vec_rModels_c_indiv[i].weight;
	}
	vec_rModels_c_indiv[0].weight=vec_rModels_c_indiv[0].weight/all_weight;
	rWeightSums_indiv.push_back(vec_rModels_c_indiv[0].weight);
	long last_item=vec_rModels_c_indiv.size()-1;
	double lci=quantile_for_CI;
	double uci=1-quantile_for_CI;
	double averaged_gamma=-199; // average only for weights between 95% CI
  	double medium_gamma=-199; // model medium gamma
  	double medium=0.5;

  	////Get the discrete model based gamma only, but not interpolation.
	//lower_model, for the case the first item is already above the CI; record the lower model as 0.  use interpolation, and the model one and model two to get the y axis intercept, to get lci gamma.
	if (rWeightSums_indiv[0]>=lci) {
		lower_model=0;
		flag_lower=1;
		double yinterp=0;

		//To Be Deleted [original way: get that position only]
		vec_lower_r_c[site]= vec_rModels_c_indiv[lower_model].r; //the Lower 95% CI gamma, not use interpolation, since another point required for that, simply use the first model as the lower 95% CI model
		averaged_gamma=(vec_rModels_c_indiv[lower_model].weight-lci)*vec_rModels_c_indiv[lower_model].r; //Get the gamma if more than the lower 95% CI
	}
  	////Added: If the first model is above the medium, record the medium_gamma
  	if (rWeightSums_indiv[0]>=medium) {
		medium_gamma=vec_rModels_c_indiv[0].r;
  	}

	//***Debug the model averaged gamma and lower CI gamma values
	//cout<<"ModelID:\t"<<"ModelAveragedGamma:\t"<<"lowerCI:\t"<< "UpperCI:\t"<<"gamma[i]:\t"<< "weight[i]:\t"<<"gamma[i-1]:\t"<<"SummedWeight"<<endl;

	for(long i=1;i<vec_rModels_c_indiv.size();i++){
		vec_rModels_c_indiv[i].weight=vec_rModels_c_indiv[i].weight/all_weight;
		double tmp_weight=vec_rModels_c_indiv[i].weight;
		if (tmp_weight == 0.0) cerr << "null weight " << site << " " << i << endl;
		rWeightSums_indiv.push_back(rWeightSums_indiv[i-1] + tmp_weight);
		////averaged_gamma+=vec_rModels_c_indiv[i].weight*vec_rModels_c_indiv[i].r;

		//Get the lower CI gamma; record the model number for the rWeightSums at the two borders of the lower 95% CI lci 0.025, lower_model
		if (rWeightSums_indiv[i-1]<=lci && rWeightSums_indiv[i]>=lci) {
			lower_model=i;
			flag_lower=1;
			vec_lower_r_c[site]= vec_rModels_c_indiv[lower_model].r;
			averaged_gamma=(rWeightSums_indiv[i]-lci)*vec_rModels_c_indiv[i].r;
			////cout<<"***Lower CI model id:\t"<< lower_model<<"\tLowerCI gamma:\t"<<lci_interpolation_rl<<"\tSummedWeight:\t"<<rWeightSums[lower_model]<<"\tGammaBefore:\t"<<vec_rModels_c_indiv[i-1].r<<"\tGammaAfter:\t"<<vec_rModels_c_indiv[i].r<<"\tAveragedGamma:\t"<<averaged_gamma<<endl<<endl;
		}
		//getting model averaged gamma between 95% CI
		if (flag_lower==1 && i>lower_model && rWeightSums_indiv[i]<=uci){
			averaged_gamma+=vec_rModels_c_indiv[i].weight*vec_rModels_c_indiv[i].r;
		}

		//Get the upper CI gamma; record the model number for the rWeightSums at the two sides of the upper 95% CI uci 0.975, upper_model
		if (flag_lower==1 && i>lower_model && rWeightSums_indiv[i-1]<=uci && rWeightSums_indiv[i]>uci) {
			upper_model=i;
			vec_upper_r_c[site]=vec_rModels_c_indiv[upper_model].r;
			averaged_gamma+=vec_rModels_c_indiv[i].r *(uci-rWeightSums_indiv[i-1]);
		}

  		////Added: medium_gamma=vec_rModels[i].r;
  		if (rWeightSums_indiv[i-1]<medium && rWeightSums_indiv[i]>=medium){
  			medium_gamma=vec_rModels_c_indiv[i].r;
  		}
		//***Debug the model averaged gamma and lower CI gamma values
		//cout<<i<<"\t"<< averaged_gamma<<"\t"<< vec_lower_r_c[site]<<"\t"<< vec_upper_r_c[site]<<"\t"<< vec_rModels_c_indiv[i].r<<"\t"<<tmp_weight<<"\t"<<vec_rModels_c_indiv[i-1].r<<"\t"<<rWeightSums[i]<<endl;
	}
	//To be deleted, model averaged gamma
	if (averaged_gamma!=-199) {
		vec_r_c[site]=averaged_gamma/(uci-lci); // model averaged gamma only takes from lci to uci, so it should divide (uci-lci).
	}
	else {
		*myout<<"Model averaged gamma: NULL"<<endl;
				vec_r_c[site]=-199;
	}
	//ZMZ 04/28/2016 use only poisson rate regional gamma as an option, by default, weighted regional gamma and recurrent gamma
	if (regional_gamma_only==0)
	{
		//Get the recurrent revised gamma and 95% gamma
		vec_r_c[site]=(1-recur_weight)*vec_r_c[site]+recur_weight*recur_r;
		vec_lower_r_c[site]=(1-recur_weight)*vec_lower_r_c[site]+recur_weight*recur_lower_r;
		vec_upper_r_c[site]=(1-recur_weight)*vec_upper_r_c[site]+recur_weight*recur_upper_r;
	}

	//Make sure the LowerCI and UpperCI is two sides of value of model average; quit the program if model averaged gamma is not between 95% CI gamma
	if (CONVERT<int>(vec_lower_r_c[site])>CONVERT<int>(vec_r_c[site]) || CONVERT<int>(vec_upper_r_c[site]) < CONVERT<int>(vec_r_c[site]) || flag_lower!=1 ||  upper_model<lower_model){
		*myout<<endl<<endl<<"Warning: Problem in getting proper values of model averaged gamma, lower and upper 95% CI gamma from all models!!!!"<<endl;
		*myout<<"Site: "<<site<<"\tGamma: "<<vec_r_c[site]<<"\tLower 95% CI gamma: "<<vec_lower_r_c[site]<<"\tUpper CI gamma: "<<vec_upper_r_c[site]<<"\tTotalModelNO: "<<last_item<<"\tLowerCIModelID:"<<lower_model<<"\tUpperCIModelID: "<<upper_model<<endl;
		*myout<<"Total number of models: "<<last_item<<" Model 1: "<<vec_rModels_c_indiv[0].r<<"\tSumWeight: "<<rWeightSums_indiv[0]<<"\tModel 2: "<<vec_rModels_c_indiv[1].r<<"\tSumWeight: "<<rWeightSums_indiv[1]<<"\tModel "<<last_item-1<<": "<<vec_rModels_c_indiv[last_item-1].r<<"\tSumWeight: "<<rWeightSums_indiv[last_item-1]<<"; Model"<<last_item<<": "<<vec_rModels_c_indiv[last_item].r<<"\tSumWeight: "<<rWeightSums_indiv[last_item]<<endl<<"Lower95% weight: "<<rWeightSums_indiv[lower_model]<<"\tUpper 95% CI weight: "<<rWeightSums_indiv[upper_model]<<endl;
		*myout<<"Converted to integers: Lower: "<<CONVERT<int>(vec_lower_r_c[site])<<"\tGamma:\t"<<CONVERT<int>(vec_r_c[site])<<"\tUpper:\t"<<CONVERT<int>(vec_upper_r_c[site])<<endl;
		*myout<< "Warning in getting the right 95% confidence interval gamma in CI_UpLow_rc!";
	}
	vec_rModels_c_indiv.clear();
	return 1;
}

/***************************************************
* Function: parse Parameters for running cMAC-PRF
* Input Parameters: polymorphism file name, divergence file name, other parameter options including
* Output: integrated input parameters
* Return Value: true or false
***************************************************/
bool cPRFCluster::parseParameter(int argc, const char* argv[]) {
  bool flag = true;
  int i;
  string temp;
	
  try {	
    if (argc==4 or argc==2) {
      temp = stringtoUpper(argv[1]); // Read in user input parameters 'temp' and make them all upper cases to be recognized below.
      if (temp=="-H") showHelpInfo();
      else {cout<<"argc number problems, 2, 4!\n"; throw 1;}
    }
    else if (argc!=5 && argc!=7 && argc!=9 && argc!=11 && argc!=13 && argc!=15 && argc!=17 && argc!=19 && argc!=21 && argc!=23 && argc!=25 && argc!=27 && argc!=29 && argc!=31 && argc!=33) {
      cout<<"argc number problems, not 5 to 33!\n";
      throw 1;			
    }		
    else {			
      //default values for parameters
      int recur_flag=0, output_flag=0, div_num_flag=0, div_cons_flag=0,code_flag=0, criterion_flag=0,ms_flag=0, synonymous_flag=0,ci_ma_flag=0,r_flag=0,Do_ci_r_flag=0,ci_method_flag=0, nuc_replace_flag=0;
      SilentRate_flag=0;
      for (i=1; i<argc; i++) {				
			temp = stringtoUpper(argv[i]);
			//option - the divergent consensus input file, required.
			if (temp=="-DC" && (i+1)<argc && div_cons_flag==0) {
			  div_cons_seqfile = argv[++i];
			  div_cons_flag++;
			}
			//option - the number of tumor samples, required.
			else if (temp=="-DN" && (i+1)<argc && div_num_flag==0) {
				tumor_num_s = argv[++i];
				div_num_flag++;
			}
			//option - the recurrent input file.
			else if (temp=="-RF" && (i+1)<argc && recur_flag==0) {
				recurrent_file = argv[++i];
			  recur_flag++;
			}
			//option - the output format, amino acid or nucleotide output; 0: amino acid output || 1: nucleotide output, default=0
			else if (temp=="-O" && (i+1)<argc && output_flag==0) {
				  int num = CONVERT<int>(argv[++i]);
				  if (num==0 || num==1) {//0: amino acid output || 1: nucleotide output, default=0			  
					  output_format_num=num;
					  output_flag++;
				  }
				  else {
					cout<<"Error in parseParameter for -O, output format should be 0 or 1, 0: amino acid output || 1: nucleotide output, default=0."<<endl;
					throw 1;
				  }
			}
			//option - user defined silent rate
			else if (temp=="-SR" && (i+1)<argc && SilentRate_flag==0) {
				SilentRate =CONVERT<double>(argv[++i]);;
				SilentRate_flag++;
			}
			//option - criteria to calculate weights, AIC/BIC/AICc/LRT
			else if (temp=="-C" && (i+1)<argc && criterion_flag==0) {
			  int num = CONVERT<int>(argv[++i]);
			  if (num>=0 && num<=3) {
				criterion_type = num;
				criterion_flag++;
			  }
			  else {
				cout<<"Error in parseParameter for -C, criteria AIC/BIC/AICc"<<endl;
				throw 1;
			  }
			}
			//Genetic code
			else if (temp=="-G" && (i+1)<argc && code_flag==0) {
			  int num = CONVERT<int>(argv[++i]);
			  if (num>0 &&  num<24) {
				genetic_code = num;
				code_flag++;
			  }
			  else {
				cout<<"Error in parseParameter for -G"<<endl;
				throw 1;
			  }
			}
			//Model selection and model averaging
			else if (temp=="-M" && (i+1)<argc &&  ms_flag==0) {
			  int num=CONVERT<int>(argv[++i]);
			  if(num==0 || num==1){
				MS_only = num;
				ms_flag++;
			  }
			  else {
				cout<<"Error for the input parameter -m, the input should be 0 or 1. {0: use both model selection and model averaging || 1: use only model selection}, default = 0}"<<endl;
				throw 1;
			  }
			}
			//Optional - Show the clustering results from synonymous sites
			else if (temp=="-S" && (i+1)<argc && synonymous_flag==0){
			  int num=CONVERT<int>(argv[++i]);
			  if(num==0 || num==1){
				Do_Synonymous_Cluster = num;
				synonymous_flag++;
			  }else{
				cout<<"Error in parseParameter for -S"<<endl;
				throw 1;
			  }
			}
			//Optional - Scale of the gene, it should be 1,3,6,9,12,15...3*n, the true sequence length is the times of the scale and the given length.
			else if (temp=="-SC" && (i+1)<argc){
			  int num=CONVERT<int>(argv[++i]);
			  Scale = num;
			  if (Scale%3!=0) { cout<< "Warning: Scale of the gene cannot be divided by 3 (codon size)."<<endl;}
			  if (num<1)
			  {
				cout<<"Error in parseParameter for -SC Scale is smaller than 1."<<endl;
				throw 1;
			  }
			}
			//option - calculation of confidence intervals for model averaging of poisson rates
			else if(temp=="-CI_M" && (i+1)<argc && ci_ma_flag==0){
			  int num=CONVERT<int>(argv[++i]);
			  if(num==0 || num==1){
					ci_ma=num;
				ci_ma_flag++;
				  }else{
					cout<<"Error in parseParameter for -CI_M."<<endl;
					throw 1;
				  }
			}
			//option - gamma estimation
			else if(temp=="-R" && (i+1)<argc && r_flag==0){
				  int num=CONVERT<int>(argv[++i]);
			  if(num==0 || num==1){
					Do_r_estimate=num;
				r_flag++;
				  }else{
					cout<<"Error in parseParameter for -R."<<endl;
					throw 1;
				  }
				}
			//option - calculation of confidence intervals of gamma
			else if(temp=="-CI_R" && (i+1)<argc && Do_ci_r_flag==0){
				  int num=CONVERT<int>(argv[++i]);	  
			  if(num==0 || num==1){
					Do_ci_r=num;
				Do_ci_r_flag++;
				  }else{
					cout<<"Error in parseParameter for -CI_R."<<endl;
					throw 1;
				  }
				}
			//option - usage of exact algorithm for estimating CIs for gamma
			else if(temp=="-EXACT" && (i+1)<argc && ci_method_flag==0){
				  int num=CONVERT<int>(argv[++i]);
			  if(num==0 || num==1){
					Do_ci_r_exact=num;
				ci_method_flag++;
				  }else{
					cout<<"Error in parseParameter for -EXACT."<<endl;
					throw 1;
				  }
				}
			//Option - regional gamma only for gamma estimation, default is weighted gamma of regional and recurrent
			//REGIONALGAMMA needs to be capitalized
			else if(temp=="-REGIONALGAMMA" && (i+1)<argc){
			  int num=CONVERT<int>(argv[++i]);
			  if(num==0 || num==1){
				regional_gamma_only=num;
				  }
			  else{
				  cout<<"Error in parseParameter for -RegionalGamma."<<endl;
				  throw 1;
				  }
			}   
			//Replace or see as gap in the ambiguous nucleotide site 
			else if(temp=="-N" && (i+1)<argc && nuc_replace_flag==0){
			  int num=CONVERT<int>(argv[++i]);
			  if(num==0 || num==1){
				Nuc_replace=num;
				nuc_replace_flag++;
			  }else{
				cout<<"Error in parseParameter for -N."<<endl;
				throw 1;
			  }
			}
			else{
				  cout<<"Error in parseParameter, not recognized parameter."<<endl;
				  throw 1;
				}
      }	//end of for argc	
    }//end of else argc
  }//end of try
  catch (...) {
    cout<<"Error in parsing input parameter(s)."<<endl;
    cout<<"Type -h for help information."<<endl;
    cout<<NAME<<", Version: "<<VERSION<<" [Last Update: "<<LASTUPDATE<<"]"<<endl;
    flag = false;		
  } 
  return flag;
}

/***************************************************
* Function: use './cMAC-PRF -h' to show help options for cMAC-PRF, especially assignments for different parameters.
* Input Parameter: none
* Output:Screen output of cMAC-PRF information and options of parameters.
* Return Value: none
***************************************************/
void cPRFCluster::showHelpInfo() {
  cout<<"***********************************************************************"<<endl;
  cout<<NAME<<", Version: "<<VERSION<<" ["<<LASTUPDATE<<"]"<<endl;	
  cout<<"Function: "<<FUNCTION<<endl;
  cout<<"Usage: "<<NAME<<" [OPTIONS]"<<endl;	
  cout<<"***********************************************************************"<<endl<<endl;
  cout<<"OPTIONS:"<<endl;	
  cout<<"  -dc\tInput file name for divergent consensus sequences [string, required]"<<endl;
  cout<<"  -dn\tInput number for divergent consensus sequences (sample number) [string, required]"<<endl;
  cout<<"  -rf\tInput replacement recurrent file [string, optional]"<<endl;
  cout<<"  -o\tChoice of output format [integer, optional], {0: amino acid level output || 1: nucleotide level output, default=0}"<<endl;
  cout<<"  -sr\tSilent Rate [preferred from MutSigCV, corrected with gene expression,replicating time,and chromosomal position.] [double, required when the number of synonymous mutation equals 0]"<<endl;
  cout<<"  -c\tCriterion used for clustering [integer, optional], {0:AIC || 1: BIC || 2:AICc || 3:LRT}, default = 0"<<endl;
  cout<<"  -g\tGenetic code used for sequences [integer, optional], {1:standard}, default = 1"<<endl;
  cout<<"    \tMore information about the genetic codes can be found at http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"<<endl;
  cout<<"  -m\tModel selection and model averaging[integer, optional], {0: use both model selection and model averaging || 1: use only model selection}, default = 0"<<endl;
  cout<<"  -ci_m\tCalculate 95% confidence intervals for results of model averaging [integer, optional], {0: NOT calculate 95% confidence intervals || 1: calculate 95% confidence intervals}, default = 0"<<endl;
  cout<<"  -s\tShow clustering results of synonymous sites from Divergent sequences [integer, optional], {0: without clustering results of synonymous and replacement sites || 1: with clustering results of synonymous and replacement sites}, default = 0"<<endl;
  cout<<"  -sc\tScale of the gene, it should be 1,3,6,9,12,15...3*n, the true sequence length is the times of the scale and the given length."<<endl;

  cout<<"  -r\tEstimate selection coefficient for each site [integer, optional], {0: NOT estimate selection coefficient || 1: estimate selection coefficient}, default=1"<<endl;
  cout<<"  -ci_r\tCalculate 95% confidence intervals for selection coefficient [integer, optional], {0: NOT calculate 95% confidence intervals || 1: calculate 95% confidence intervals}, default = 1"<<endl;
  cout<<"  -exact\tAlgorithm for calculating 95% confidence intervals for selection coefficient [integer, optional], {0: use stochastic algorithm || 1: use exact algorithm}, default = 0"<<endl;
  //ZMZ 04/28/2016  added option - regional gamma only
  cout<<"  -RegionalGamma\tRegional gamma only [integer, optional], {0: use weighted regional and recurrent gamma || 1: use regional gamma only}, default = 0"<<endl;

  cout<<"  -h\tShow help information"<<endl;	
  cout<<endl;
	
  cout<<"COPYRIGHT & LICENSE:"<<endl;
  cout<<NAME<<" is distributed as open-source software and licensed under the GNU General Public License "<<endl;
  cout<<"(Version 3; http://www.gnu.org/licenses/gpl.txt), in the hope that it will be useful, but "<<endl;
  cout<<"WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A "<<endl;
  cout<<"PARTICULAR PURPOSE. See the GNU General Public License for more details."<<endl;
  cout<<endl;
  
  cout<<"REFERENCE:"<<endl;	
  cout<<REFERENCE<<endl;
  cout<<endl;

  cout<<"SEE ALSO:"<<endl;
  cout<<"For more information, please see <http://www.yale.edu/townsend/software.html>."<<endl;
  cout<<endl;
  
  cout<<"CONTACT:"<<endl;
  cout<<"Please send suggestions or bugs to Ziming Zhao (ziming.gt@gmail.com) and Jeffrey Townsend (Jeffrey.Townsend@Yale.edu)."<<endl;
  cout<<endl;
}
/***************************************************
Subfunction - open and read the lookup table, to extract the values for k - the count of recurrent mutations for a particular site.
Input: the lookup table file name, the struct LambdaCI;
Output: vector<CI_recurrent> RecurrentLookup LambdaCI
***************************************************/
int cPRFCluster::LambdaCILookupTable(string input_f_name){
	int k=0;
	double tmp_l=0.0, tmp_u=0.0;
   ifstream myfileFn2(input_f_name.c_str());
   if (!myfileFn2) throw "Error in opening LambdaCILookupTable...\n";
   //Read the file; remove the empty lines at the end of the file
   string str;
   getline(myfileFn2,str);
   	while ( myfileFn2.good()) {
   	   getline(myfileFn2,str);
   	  // cout<<"Each line: "<<str<<endl;
   	  int position1 = str.find("	");
   	     if (position1!=std::string::npos)
   	     {
   	    	 //cout<<"The position1: "<<position1<<endl;
   	     }
   	     else {cout<<"Error! Failed to find the position1 in the lookup table!"<<endl;}
   	     int position2=str.find("	", position1+2);
		  string e1 = str.substr (0,position1);
		string e2 = str.substr (position1+1,position2-position1-1);
		string e3 = str.substr (position2+1,str.length()-position2-1);
		//cout<<e1<<"\t**"<<e2<<"$$\t***"<<e3<<endl;
		k=CONVERT<int>(e1);
		tmp_l= CONVERT<double>(e2);
		tmp_u=CONVERT<double>(e3);
	   CIRecurrentLookup tmp_ci(k,tmp_l,tmp_u);
	   LambdaCIs.push_back(tmp_ci);
	   //cout<<k<<"**\t"<<tmp_l<<"$$\t***"<<tmp_u<<endl;
	   k=0; tmp_l=tmp_u=0.0;
   	}
   	myfileFn2.close();
   	//cout<<"***The size of the lookup LambdaCI: "<<LambdaCIs.size()<<endl;
    //cout<<"***Test Lambda lookup table:"<<endl;
    //cout<<LambdaCIs[0].lowerCIlambda<<endl;
    //cout<<LambdaCIs[0].count<<endl;
 return 1;
}

/***************************************************
* Function: Estimate parameter selection intensity r for recurrent site for cancer, and CI for lowerCILambda and upperCILambda
* Input Parameter: the number of species, the mutation rate, the recurrent num, and site
* Output: gamma for recurrent site
* Return Value: gamma for recurrent site
***************************************************/
int cPRFCluster::recurrent_SiteGamma(int tumor_num, double ucr, int RecurrentNum, int Site, ostringstream* myout){
	double j=1.0;
	//Calculate the parameter r for recurrent sits
	double tmp=(ucr/RecurrentNum)*tumor_num;
	double df,f;
	double rtn,dx;

	double min_dx=IR_H;
	double optimal_r=IR_H;

	bool flag_root=false;
	for(rtn=IR_L;rtn<=IR_H;rtn+=0.5){
		if(flag_root==true){
			break;
		}
		for(j=1; j<=JMAX; j++){
			double tmp0=-2*rtn; //rtn is gamma
			f=(1-exp(tmp0))/(2*rtn)-tmp;//f(r0) after trapezoidal rule numeric integration of x
			df=(2*rtn*exp(tmp0)-(1-exp(tmp0)))/(2*rtn*rtn);//f(r0)'after trapezoidal rule numeric integration of x
			dx=f/df; //f(r0)/f(r0)'
			rtn -= dx; // iteration for gamma in Newton's method, r1=r0-f(r0)/f(r0)'
			if(fabs(dx)<fabs(min_dx)){
				min_dx=dx;
				optimal_r=rtn;
			}
			//criteria to quit the iterations for gamma calculation.
			if(fabs(dx)<ER){
				//vec_r_c[Site]=rtn;
				vec_r_c_r[Site]=rtn;
				flag_root=true;
				break;
			}
		}
	}
	if(flag_root==false and min_dx<=MinDx){
		vec_r_c_r[Site]=optimal_r;
	}
	if(flag_root==false and min_dx>MinDx){
		//vec_r_c[Site]=-199;
		vec_r_c_r[Site]=-199;
	}
	*myout <<"Calculate the parameter r for recurrent site "<<Site<<" with Recurrent number "<<RecurrentNum<<" Gamma: "<<vec_r_c_r[Site]<<endl;
	return 1;
}

/***************************************************
* Function: Estimate parameter selection intensity r for recurrent site for cancer, and CI for lowerCILambda and upperCILambda
* Input Parameter: the number of species, the mutation rate, the recurrent num, and site
* Output: gamma for recurrent site
* Return Value: gamma for recurrent site
***************************************************/
int cPRFCluster::recurrent_SiteGammaCI(int tumor_num, double ucr, double CI, int Site, int upORlow, ostringstream* myout){
	double j=1.0;
	//Calculate the parameter r for recurrent sits
	double tmp=(ucr/CI)*tumor_num;
	double df,f;
	double rtn,dx;
	double tmp_gamma;
	double min_dx=IR_H;
	double optimal_r=IR_H;
	bool flag_root=false;
	for(rtn=IR_L;rtn<=IR_H;rtn+=0.5){
		if(flag_root==true){
			break;
		}
		for(j=1; j<=JMAX; j++){
			double tmp0=-2*rtn; //rtn is gamma
			f=(1-exp(tmp0))/(2*rtn)-tmp;//f(r0) after trapezoidal rule numeric integration of x
			df=(2*rtn*exp(tmp0)-(1-exp(tmp0)))/(2*rtn*rtn);//f(r0)'after trapezoidal rule numeric integration of x
			dx=f/df; //f(r0)/f(r0)'
			rtn -= dx; // iteration for gamma in Newton's method, r1=r0-f(r0)/f(r0)'
			if(fabs(dx)<fabs(min_dx)){
				min_dx=dx;
				optimal_r=rtn;
			}
			//criteria to quit the iterations for gamma calculation.
			if(fabs(dx)<ER){
				tmp_gamma=rtn;
				flag_root=true;
				break;
			}
		}
	}
	if(flag_root==false and min_dx<=MinDx){
		tmp_gamma=optimal_r;
	}
	if(flag_root==false and min_dx>MinDx){
		tmp_gamma=-199;
	}
	if (upORlow==0) {
		//vec_lower_r_c[Site]=tmp_gamma;
		vec_lower_r_c_r[Site]=tmp_gamma;
		*myout<<"Lower recurrent gamma: "<<vec_lower_r_c_r[Site]<<endl;
	}
	else {
		//vec_upper_r_c[Site]=tmp_gamma;
		vec_upper_r_c_r[Site]=tmp_gamma;
		*myout<<"Upper recurrent gamma: "<<vec_upper_r_c_r[Site]<<endl;
	}
	return 1;
}
/***************************************************
* Function: Get a random model by choosing a value from 0 to 1, and select the clustering model for p with the weight above the random value
* Input Parameter: SiteModels; site
* Output:model_num
* Return Value: model_num
***************************************************/
long cPRFCluster::RandomModel_Num(struct SiteModels *p,long site){
  long model_num=-1;
  //randomize();
  double rand_tmp=rand();
  //RAND_MAX is from system (2147483647)
  double rand_num=rand_tmp/RAND_MAX;
  if(rand_num==0.0){
      return 0;
  }
  //double rand_tmp=rand()%(N_Random*N_Random);
  //double rand_num=rand_tmp/(double)(N_Random*N_Random);
  int pss = p[site].sms.size();
  double sum=0.0;
  CI *cip = &(p[site].sms[0]);
  for(long i=0;i<pss;i++,cip++){
    double pwt = cip->weight;
    if(rand_num>sum && rand_num <= pwt+sum){
      model_num=i;
      break;
    }else{
      sum=pwt+sum;
    }
  }
  return(model_num);
}

//BUG: The model averaged gamma should not be the gamma corresponding to the model averaged probability; it should be gamma averaged for all gamma calculated from all probabilities.
/***************************************************
* Function: Estimate parameter selection intensity r for each site for cancer
* Input Parameter: the number of species, the sequence length
* Output: gamma for each site
* Return Value: gamma for each site
***************************************************/
int cPRFCluster::rc_SitePRF(int tumor_num, double ucr, long N){
  double j=1.0;
  //Estimate r for each site one by one
  for(long i=0; i<N; i++){
	  double min_dx=IR_H;
	  double optimal_r=IR_H;

    if(vec_MA_rate_dr[i]==0){
      //Under infinite negative selection
      vec_r_c[i]=0; //ZMZ 07/18/2016 
      continue;
    }else if(ucr==0){
      vec_r_c[i]=299;
      continue;
    }

    //Calculate the parameter r
    double tmp=(ucr/vec_MA_rate_dr[i])*tumor_num;
    double df,f;
    double rtn,dx;

    bool flag_root=false;
    for(rtn=IR_L;rtn<=IR_H;rtn+=0.5){
    	if(flag_root==true){
    		break;
    	}
    	for(j=1; j<=JMAX; j++){
    		double tmp0=-2*rtn; //rtn is gamma
    		f=(1-exp(tmp0))/(2*rtn)-tmp;//f(r0) after trapezoidal rule numeric integration of x
    		df=(2*rtn*exp(tmp0)-(1-exp(tmp0)))/(2*rtn*rtn);//f(r0)'after trapezoidal rule numeric integration of x
    		dx=f/df; //f(r0)/f(r0)'
    		rtn -= dx; // iteration for gamma in Newton's method, r1=r0-f(r0)/f(r0)'
    		if(fabs(dx)<fabs(min_dx)){
    			min_dx=dx;
    			optimal_r=rtn;
    		}
    		//criteria to quit the iterations for gamma calculation.
    		if(fabs(dx)<ER){
    			vec_r_c[i]=rtn;
    			flag_root=true;
    			break;
    		}
    	}
    }
    if(flag_root==false and min_dx>MinDx){
    	cout<<"rs_SitePRF: Site:\t"<<i<<"\tGamma NULL. min_dx:\t"<<min_dx<<"\tOptimalGamma:\t"<<optimal_r<<"\tRate:\t"<<vec_MA_rate_dr[i]<<"\tucr:\t"<<ucr<<"\tTumorNumber:\t"<<tumor_num<<endl;
    	vec_r_c[i]=-199;
    }
    if(flag_root==false and min_dx<=MinDx){
    	vec_r_c[i]=optimal_r;
    }
    //cout<<"Gamma: "<<vec_r_c[i]<<"\toptimal_r: "<<optimal_r<<"\tmin_dx: "<<min_dx<<endl;
  }
  return 1;
}