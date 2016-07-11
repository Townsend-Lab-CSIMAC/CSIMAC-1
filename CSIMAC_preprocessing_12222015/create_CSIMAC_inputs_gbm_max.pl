#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  create_cMACPRF_inputs.pl
#
#        USAGE:  perl create_cMACPRF_inputs.pl Combined_LUSC 2 2 mutationsTN_23_Lung_Squamous_Cell_Carcinoma.maf mutationsTN_21_TCGA_LUSC.maf;
#                perl create_cMACPRF_inputs.pl Gilead_LUSC 2 1 mutationsTN_23_Lung_Squamous_Cell_Carcinoma.maf;
#                perl create_cMACPRF_inputs.pl Gilead_LUAD 2 1 mutationsTN_22_Lung_adenocarcinoma.maf 
#                perl create_cMACPRF_inputs.pl Combined_LUAD 2 2 mutationsTN_22_Lung_adenocarcinoma.maf mutationsTN_18_TCGA_LUAD.maf
#                perl create_cMACPRF_inputs.pl Gilead_Melanoma 2 1 mutationsTN_1_Melanoma_exome_sequencing.maf;
#perl create_cMACPRF_inputs.pl luad_both_trim 2 1 luad_combined_2013-12-02_trim_SNP.maf
#perl create_cMACPRF_inputs.pl Meningiomas 2 1 CancerDB_2014-02-17_meningioma_maf_yg_format.maf.tsv;
#  DESCRIPTION: Creates input for cMACPRF which is a string of * with R, S and B; R-replacement, S-Synonymous, B-both (mutations) 
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  Fixed 11/05/2013 The pbs file, the recurrent file name - removing the directory. 
#        NOTES:  ---
#       AUTHOR:   (Zi-Ming Zhao), <ziming.gt@gmail.com>
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  03/20/2013 10:53 AM EDT
#      UPDATED:  11/05/2013
#===============================================================================
#Attention before each run: change the gene list for the specific tumor type
#Attention before pbs runs: need LookupTable_cMACPRF_CI_Recurrent_v9.dat in the folder
#Attention: For Melanoma, $file_number2=0 causes problem, thus both file numbers are one larger than it should be.
	#Get the overall genomic mutation rate and the smallest mutation rate
	#my $inf_overall=$infolder.$TumorType.".mutsigcv_overall_rates.txt";
	
# Oct 06, 2015 Ovary Cancer pre-run Ziming
# Command line: perl create_CSIMAC_inputs_gbm_max.pl Max_GBM_Data_Named 2 1 Max_GBM_Data_Named.maf 
# To revise each time: directory $dir,$indir, $outdir 
# To make sure have all input files: UCSC_refseq_refFlat_hg19; maf, mutsigcv two files: gene_rates.txt, overall_rates.txt; pbs

# Command line: perl create_CSIMAC_inputs_gbm_max.pl Max_GBM_Data_Named 2 1 Max_GBM_Data_Named.maf >gbm_screenOutput.txt &


#use strict;
use warnings;
use Carp;
use Data::Dumper;
#use GFF;
#use RefFlat;
use MAF_ziming;
use Math::Cephes qw(:dists);
use Tie::SortHash;
use List::MoreUtils qw(uniq);
#use Tie::IxHash;

# Change the scale based on the gene length
my $start_time = time();

#my $dir="/Users/zimingzhao/Desktop/";
my $dir = "/Users/zimingzhao/Desktop/XinruRen/Pre-processing_gbm_max/";
my $indir=$dir;
#my $indir = $dir."/Pre-processing/";
#my $outdir = $dir."cMAC-PRF/cMACPRF_inputs/";
my $outdir = $dir."Pre-processing_Output_max_v2/";

# UCSC refFlat file downloaded from UCSC Table Browser
print "\n\n......Reading the file with cds and exon information used - UCSC_refseq_refFlat_hg19......\n";
my $USAGE2 = "\nPlease provide file with cds and exon information used in UCSC_refseq_refFlat_hg19\n";
#my $path_to_genomeinfo = '/Users/zimingzhao/Desktop/GenomicInfo/';
my $path_to_genomeinfo = $dir;
my $refFlat=$path_to_genomeinfo."UCSC_refseq_refFlat_hg19";
my @RefSeq=@{MAF::parseREFSEQ($refFlat)};

# Alternative using RefFlat package
#my $refFlat_obj = RefFlat->new($refFlat);
#print "RefFlat file name: ", $refFlat_obj->get_refFlat_filename, "\n";
#my @RefFlatData = @{$refFlat_obj->get_refFlat_data};
 
print "\n\n......Reading tumor type, the numbers of file names (Yale; TCGA), and MAF files (Yale; TCGA) as inputs......\n";
my $USAGE1 = "\nPlease, provide tumor type, the numbers of file names (Yale; TCGA), and MAF files (Yale) as inputs\n";
my $USAGE3 = "\nPlease, provide tumor type, the numbers of file names (Yale; TCGA), and MAF files (TCGA) as inputs\n";
my $tumor_type= $ARGV[0] or croak $USAGE1; # as part of the output
my $file_number1= $ARGV[1] or croak "\nPlease enter file number 1\n"; # The number of input files
printf "Tumor type: %s; File number 1: %d;\n",$tumor_type,$file_number1;
my $file_number2= $ARGV[2] or croak "\nPlease enter file number 2\n"; # The number of input files

my @mutations;
my @allgeneids;
printf "Tumor type: %s; File number 1: %d; File number 2: %d\n",$tumor_type,$file_number1,$file_number2;
print "\n\n......opening MAF files with mutations......\n";

for (my $ii=3;$ii<$file_number1+2;$ii++)
	{
    	# my $single_mut_file = $indir.$ARGV[$ii] or croak $USAGE1;
    	my $single_mut_file = $indir."Max_GBM_Data_Named.maf" or croak $USAGE1;
    	my ($loc_tmp_mutations, $loc_geneids) = MAF::parse_MAF3333($single_mut_file) ; # Yale
    	my @tmp_mutations = @{$loc_tmp_mutations};
    	my @geneids = @{$loc_geneids};
    	push @mutations, @tmp_mutations; 
    	push @allgeneids, @geneids;
	}

for (my $jj=$file_number1+2;$jj<$file_number2+$file_number1+1;$jj++)
	{
    	my $single_mut_file = $indir.$ARGV[$jj] or croak $USAGE3;
    	my ($loc_tmp_mutations, $loc_geneids) = MAF::parse_MAF3333($single_mut_file) ; # TCGA (parse_MAF33 is TCGA)
    	my @tmp_mutations = @{$loc_tmp_mutations};
    	my @geneids = @{$loc_geneids}; 
    	push @mutations, @tmp_mutations;
        push @allgeneids, @geneids;
	}
	
print scalar(@mutations), " mutations in the list.\n";
print scalar(@allgeneids), " gene IDs in the list.\n";

#Gene list to run in CSIMAC
#my @GeneIDs=qw(TP53);
my @GeneIDs=MAF::RemoveRecurrent(@allgeneids);
print scalar(@allgeneids), " gene IDs after removing redundant genes.\n";

my $CSIMAC_simplequeue = "cd ~/scratch/CSIMAC/TumorType; /home2/ad938/bin/CSI-MAC -DC TumorcMACPRFMutInput.txt -DN TumorSampleNumber -S 1 -SC Scale -SR BackgroundMutationRate -RF TumorcMACPRFRecurrentInput.txt >GeneID.TumorType_louise.CSIMAC.output.txt;\n";

my $qsub_all="";
my $cMACPRF_all="";
my $commands = "";
my @NotFoundRates;
my @ZeroRates;
my @FoundRates;

## GAH Edit
my $summary="";
my @NotSelected="";
my @Damaged = "";
## GAH Edit

#Get tumor sample number
my $TumorNumber=MAF::TumorNumberCount(\@mutations);
$summary.="Gene ID\tGeneLength\tReplacements\tNS_Recurrent\tSyn_Recurrent\n";
        
foreach my $geneid (@GeneIDs)
	{
        my $incommand = $CSIMAC_simplequeue;
		if (Create_cMACPRFinput($indir,$outdir,$TumorNumber,\@RefSeq,\@mutations,$tumor_type,$geneid,$incommand) eq "ERROR"){
		printf "\n\n\n*******Warning: Gene %s was not found any isoforms in the RefFlat!\n",$geneid;
		next;
		}
		
		if (Create_cMACPRFinput($indir,$outdir,$TumorNumber,\@RefSeq,\@mutations,$tumor_type,$geneid,$incommand) eq "SKIP") {
		next;
		}
		
		my($qsub,$cMACPRF,$MutSigCVstatus,$final_length,$replacements,$nrecurrent,$srecurrent,$notrelevant,$damaged,$command)=Create_cMACPRFinput($indir,$outdir,$TumorNumber,\@RefSeq,\@mutations,$tumor_type,$geneid,$incommand);
		$qsub_all.=$qsub;
		$cMACPRF_all.=$cMACPRF;
        $commands.=$command;
        
        #GAH Edit
        $summary.=$geneid."\t".$final_length."\t".$replacements."\t".$nrecurrent."\t".$srecurrent."\n";
        
        if ($notrelevant == 1) {
            push @NotSelected,$geneid;
            print "\nNot relevant.\n";
            next;
        }
        if ($damaged == 1) {
            push @Damaged,$geneid;
            print "\nGene damaged!\n";
            next;
        }
        if (!defined my($MutSigCVstatus)) { print "MutSigCV status undefined!"; next;}
        ##GAH Edit
        if (my ($MutSigCVstatus) eq "NotFound") { push @NotFoundRates,$geneid;}
		elsif (my ($MutSigCVstatus) eq "Zero") { push @ZeroRates,$geneid;}
		elsif (my ($MutSigCVstatus) eq "Regular") { push @FoundRates,$geneid;}
		else { croak "Error in finding right MutSigCV status for the gene ",$geneid,"!\n";}
	}

printf "\n\n***********\nMutSigCV gives %d regular background mutation rates, %d zero rates, %d genes not covered by MutSig!*********\n\n",scalar(@FoundRates),scalar(@ZeroRates),scalar(@NotFoundRates);
if (scalar(@FoundRates)<20 && scalar(@FoundRates)>0) { printf "\nThe list of genes which had mutsigcv rates are: \n "; print join "\t",@FoundRates; print "\n"; }
if (scalar(@ZeroRates)<20 && scalar(@ZeroRates)>0) { printf "\nThe list of genes with zero mutsigcv rates are: \n "; print join "\t",@ZeroRates;  print "\n"; }
if (scalar(@NotFoundRates)<20 && scalar(@NotFoundRates)>0) { printf "\nThe list of genes not covered by mutsigcv are: \n "; print join "\t",@NotFoundRates; print "\n";  }

print "There are: ", scalar(@allgeneids), " gene IDs in the gene list.";

###GAHEdit
if (scalar(@NotSelected)<20 && scalar(@NotSelected)>0) { printf "\nThe genes that were not selected are:"; print join "\t", @NotSelected;}
    #print "\n"; }
if (scalar(@Damaged)<20 && scalar(@Damaged)>0) { printf "\nThe damaged genes are:"; print join "\t", @Damaged; print "\n\n"; }
###GAHEdit

MAF::WriteFile($outdir,"qsub_pbs_cMACPRFcommands_".$tumor_type.".txt",$qsub_all);
MAF::WriteFile($outdir,"cMACPRFcommands_".$tumor_type.".txt",$cMACPRF_all);
MAF::WriteFile($outdir, "summary.txt", $summary);
MAF::WriteFile($outdir, "simplequeue.txt", $commands);

my $end = time();
my $elapsed_time = $end - $start_time;
print "Elapsed Timed  $elapsed_time\"\n";
exit;

###############################################################################
## Subroutines
###############################################################################
sub Create_cMACPRFinput
{
	my ($infolder,$outfolder,$TumorNumber,$refseq_ref,$mutation_ref,$tumor_type,$geneid, $command) = @_;
    
	#Extracting all isoforms' coordinates; if empty, return 0, go out

	my ($longest_isoform_length,$strand,$accession_ref,$ref_coordinates_for,$ref_isoform_length_for) = MAF::get_all_isoforms_refflat($geneid,$refseq_ref);
	if (MAF::get_all_isoforms_refflat($geneid,$refseq_ref)==0) {return "ERROR";}
	
	my %isoform_length_for=%{$ref_isoform_length_for};
	
	#Extracting all mutations for the given geneid from MAF files; convert genomic positions to CDS positions.
	my ($longest_isoform_syn,$max_isoform_syn,$ref_trans_tumor_syn,$longest_isoform_nonsyn,$max_isoform_nonsyn,$ref_trans_tumor_nonsyn,$longest_isoform_damaging,$max_isoform_damaging,$ref_trans_tumor_damaging)=MAF::ExtractMutPos($strand,$ref_coordinates_for,$ref_isoform_length_for,$geneid,$mutation_ref);
	my %trans_tumor_nonsyn=%{$ref_trans_tumor_nonsyn};
	my %trans_tumor_syn=%{$ref_trans_tumor_syn};	
	my %trans_tumor_damaging=%{$ref_trans_tumor_damaging};
    
	#Select the longest isoform with most missense mutations, etc....
	my $final_isoform=$max_isoform_nonsyn;
	my $final_length=$isoform_length_for{$final_isoform};
	
    my @final_tumor_nonsyn_pos= @{$trans_tumor_nonsyn{$final_isoform}};
	my @final_tumor_syn_pos= @{$trans_tumor_syn{$final_isoform}};
	my @final_tumor_damaging_pos= @{$trans_tumor_damaging{$final_isoform}};
	
	#Get the redundant positions and counts, as the cMACPRF inputs		
	my ($ref_redun_silent)=MAF::KeepRecurrent(@final_tumor_syn_pos);
    my ($ref_redun_missense)=MAF::KeepRecurrent(@final_tumor_nonsyn_pos);
	my ($ref_redun_damaging)=MAF::KeepRecurrent(@final_tumor_damaging_pos);
    
    my @nonredun_final_tumor_syn_pos= MAF::RemoveRecurrent(@final_tumor_syn_pos);
	my @nonredun_final_tumor_nonsyn_pos= MAF::RemoveRecurrent(@final_tumor_nonsyn_pos);
	my @nonredun_final_tumor_damaging_pos= MAF::RemoveRecurrent(@final_tumor_damaging_pos);
	#count the number of replacement mutations and keep only total replacement mutations 
	my $count_ns=scalar(@final_tumor_nonsyn_pos);
	my $count_ns_nonrecur=scalar(@nonredun_final_tumor_nonsyn_pos);
	if ($count_ns<2) { 
      printf "Skip gene %s with total replacement: %d\n",$geneid,$count_ns;
      return "SKIP";
	}
	
		
	#Scale the sequence mutation positions and redundant sites based on gene length
  my $scale=1;
  if ($final_length>900 and $final_length<=2700) {$scale=3;}
  elsif ($final_length>2700 and $final_length<=5400) {$scale=6;}
  elsif ($final_length>5400 and $final_length<=8100) {$scale=9;}
  elsif ($final_length>8100 and $final_length<=10800) {$scale=12;}
  elsif ($final_length>10800 and $final_length<=13500) {$scale=15;}
  elsif ($final_length>13500 and $final_length<=16200) {$scale=18;}
  elsif ($final_length>16200 and $final_length<=18900) {$scale=21;}
  elsif ($final_length>18900 and $final_length<=21600) {$scale=24;}
  elsif ($final_length>21600 and $final_length<=24300) {$scale=27;}
  elsif ($final_length>24300 and $final_length<=27000) {$scale=30;}
  elsif ($final_length>60000 and $final_length<=29700) {$scale=33;}
  elsif ($final_length>29700 and $final_length<=32400) {$scale=36;}
  elsif ($final_length>32400 and $final_length<=35100) {$scale=39;}
  elsif ($final_length>35100 and $final_length<=37800) {$scale=42;}
  elsif ($final_length>37800 and $final_length<=40500) {$scale=45;}
  elsif ($final_length>40500 and $final_length<=43200) {$scale=48;}
  elsif ($final_length>43200 and $final_length<=45900) {$scale=51;}
  elsif ($final_length>45900 and $final_length<=48600) {$scale=54;}
  elsif ($final_length>48600 and $final_length<=51300) {$scale=57;}
  elsif ($final_length>51300 and $final_length<=54000) {$scale=60;}
  elsif ($final_length>54000 and $final_length<=56700) {$scale=63;}
  elsif ($final_length>56700 and $final_length<=59400) {$scale=66;}
  elsif ($final_length>59400 and $final_length<=62100) {$scale=69;}
  elsif ($final_length>62100 and $final_length<=64800) {$scale=72;}
  elsif ($final_length>64800) {croak "The gene $geneid is too long $final_length, put a relevant scale number!\n";}


	#Creating RS sequences for cMACPRF input and recurrent input files       
  #my @tumor_scaled=CreateRS_ScaledSeq($scale,$geneid."_tumor_scaled".$scale,$final_isoform,$final_length,\@nonredun_final_tumor_syn_pos,\@nonredun_final_tumor_nonsyn_pos);
    
    #GAH Edit Check the size (SCALAR) of list of @nonredun_final_tumor_nonsyn_pos
    my @tumor_scaled=CreateRSD_ScaledSeq($scale,$geneid."_tumor_scaled".$scale,$final_isoform,$final_length,\@nonredun_final_tumor_syn_pos,\@nonredun_final_tumor_nonsyn_pos,\@nonredun_final_tumor_damaging_pos);

	my $mut_file_name=$geneid."_mut_scaled".$scale."_".$tumor_type;
	my $final_mut_file_name=PrintSeq($outfolder,\@tumor_scaled,$mut_file_name,$final_length,$final_isoform);
   
    #GAH Edit
    #Check the number of lines in this file to get RECURRENT mutations
    #see other perl code
	my $recur_file_name_nonsyn=$geneid."_nonsyn_scaled".$scale."_".$tumor_type;
	my ($final_nonsyn_redun_file_name, $nonsyn_recur_count) = PrintSeqRedun($scale,$ref_redun_missense,$recur_file_name_nonsyn,$final_length,$final_isoform);
    
    my $recur_file_name_syn=$geneid."_syn_scaled".$scale."_".$tumor_type;
	my ($final_syn_redun_file_name, $synon_recur_count)= PrintSeqRedun($scale,$ref_redun_silent,$recur_file_name_syn,$final_length,$final_isoform);

	my $recur_file_name_damaging=$geneid."_damaging_scaled".$scale."_".$tumor_type;
	my ($final_damaging_redun_file_name, $damag_recur_count) = PrintSeqRedun($scale,$ref_redun_damaging,$recur_file_name_damaging,$final_length,$final_isoform);
    
	#Get gene specific background mutation rate from MutSigCV output
	my ($MutRate,$MutsigCVstatus)=GetMutRate($infolder,$geneid,$tumor_type);
	
	#Create pbs batch files and qsub pbs running commands
    my ($qsub,$cMACPRF)=MAF::Create_pbs_batch($infolder,$outfolder,$scale,$MutRate,$TumorNumber,$tumor_type,$geneid,$final_mut_file_name,$final_nonsyn_redun_file_name);
    
    #my $command = "cd ~/scratch/CSIMAC/TumorType; /home/gah32/bin/CSI-MAC -DC TumorcMACPRFMutInput.txt -DN TumorSampleNumber -S 1 -SC Scale -SR BackgroundMutationRate -RF TumorcMACPRFRecurrentInput.txt >GeneID.TumorType_bulldogn.cMACPRF.output.txt;\n";
    
    $command =~ s/GeneID/$geneid/g;
    $command =~ s/TumorcMACPRFMutInput.txt/$final_mut_file_name/g;
    $command =~ s/TumorcMACPRFRecurrentInput.txt/$final_nonsyn_redun_file_name/g;
    $command =~ s/TumorType/$tumor_type/g;
    $command =~ s/TumorSampleNumber/$TumorNumber/g;
    $command =~ s/BackgroundMutationRate/$MutRate/g;
    $command =~ s/Scale/$scale/g;
    
    $replacements = scalar(@nonredun_final_tumor_nonsyn_pos);
    
    #GAH Edit
    my $nrecurrent = $nonsyn_recur_count;
    my $srecurrent = $synon_recur_count;
    my $notrelevant = 0;
    my $damaged = 0;

    if ($damag_recur_count > 0) {
        $damaged = 1;
    }
    ## If @final_tumor_nonsyn_pos (replacement) has less than 2 AND if location for missense recurrent mutations is empty, set return values to undef
    
    if (scalar @final_tumor_nonsyn_pos < 2 and $recurrent == 0)  {
        $notrelevant = 1;
    }
    ######### GAH Edit
	return ($qsub,$cMACPRF,$MutsigCVstatus,$final_length,$replacements,$nrecurrent,$srecurrent,$notrelevant, $damaged, $command);
}
###############################################################################



###############################################################################
#From MutSigCV output, Get the gene specific mutation for each gene; if the gene does not exist, use the genomic average; if the rate is 0, use the smallest genomic rate 
sub GetMutRate
{
	my ($infolder,$geneid,$TumorType)=@_;

	#Get the overall genomic mutation rate and the smallest mutation rate
	my $inf_overall=$infolder.$TumorType.".mutsigcv_overall_rates.txt";
	my $smallestMutRate=-1;
	my $GenomicAverageMutRate=-1;
	my $status="Regular";
	
	my $overall=MAF::openf($inf_overall);
	my $loc1=index($overall,"Silent+noncoding rate:	");
	my $loc11=index($overall,"\n",$loc1+1);
	my $loc2=index($overall,"Min non-zero x/X:	");
	my $loc22=index($overall,"\n",$loc2+1);
	
	if ($loc1!=-1 and $loc11!=-1 and $loc2!=-1 and $loc22!=-1)
		{
		 $smallestMutRate=substr($overall,$loc2+18,$loc22-$loc2-18);
		 $GenomicAverageMutRate=substr($overall,$loc1+23,$loc11-$loc1-23);
		}
	else { croak "Error in getting the smallest and genomic average mutation rate in subroutine GetMutRate!\n";}		
	
	#parse the mutation rate for genes
	my $inf_gene=$infolder.$TumorType.".mutsigcv_gene_rates.txt";		
	my @GeneRate=@{ MAF::MutsigCVparse($inf_gene)};		
	my $rate=-1;
	FEATURE:
	foreach my $datum (@GeneRate) {		
		my $gene= $datum->{gene};
		next if ($gene !~ $geneid);  
		$rate= $datum->{r_x_X};		
		}
	if ($rate==0) {
		$rate=$smallestMutRate; 
		$status="Zero"; 
		printf "Gene %s has zero background mutation!\n",$geneid;
		}
	elsif ($rate==-1) {
		$rate=$GenomicAverageMutRate; 
		$status="NotFound"; 
		printf "Gene %s is not found in MutsigCV file!\n",$geneid;
		}	
	return ($rate,$status);
}
###############################################################################



###############################################################################
sub PrintSeqRedun{
	my ($scale,$ref_new,$GeneID,$length,$ID)=@_;
	my %tumor_redun=%$ref_new;
	my %tumor_redun_scaled;
	my @new_pos;
	my $i=0;

	#scale the recurrent sites, 	
	#Attention: after scaling, new recurrent sites might happen
	foreach my $pos (sort keys %tumor_redun) {
	next if ($pos eq "");
	$i++;
	my $scaled_pos=int($pos/$scale);
	my $count=$tumor_redun{$pos};
	#save the first one
	#printf "Positions before scaling is %d\n",$pos;
	if ($i==1)
		{	
			$tumor_redun_scaled{$scaled_pos}=$count;
			push @new_pos,$scaled_pos;	
		}
	#check each nonfirst one, add on the same pos if existing already	
	elsif ( $i>1 and grep { $_->[0] eq $scaled_pos } @new_pos ) 
	
		{
			printf "******Position %d: Redundant after scaling!!!\n",$scaled_pos;
			$count=$count+$tumor_redun_scaled{$scaled_pos};
			$tumor_redun_scaled{$scaled_pos}=$count;	
    
		}
	#save the list if the scaled position is not recurrent	
	else
		{
			$tumor_redun_scaled{$scaled_pos}=$count;
			push @new_pos,$scaled_pos;		
		}

	
	}
    printf "******Before scaled redundants:\n";
    print Dumper (\%tumor_redun);
    printf "******After scaled redundants:\n";
    print Dumper (\%tumor_redun_scaled);
    
	my $output_file_name=$outdir.$GeneID."_".$length."_".$ID."_recurrent_".$i.".txt";
	my $return_output_file_name=$GeneID."_".$length."_".$ID."_recurrent_".$i.".txt";
	open( OUT1, ">".$output_file_name) or die $!;
	my $pres1 = select(OUT1);	
	
    #Return the size of the %tumor_redun_scaled hash
	foreach my $pos (sort keys %tumor_redun_scaled)
	{
	 print $pos, " => ",$tumor_redun_scaled{$pos},"\n";
	}
	
    #print Dumper (\%tumor_redun_scaled);
	select $pres1;
	close(OUT1) or die $!;
	return ($return_output_file_name, scalar keys %tumor_redun_scaled);
}
###############################################################################

###############################################################################
# make sure the position is exactly matching, need to subtract 1 for starting at 0
# $Replacement[$j]=$Replacement[$j]-1; 
sub CreateRSD_ScaledSeq
{
	my ($scale,$geneID,$isoform,$length1,$ref_Synonymous,$ref_Replacement,$ref_damaging)=@_;

	my @new2;

	my $length=int($length1/$scale);
	for (my $i=0;$i<$length;$i++)
		{
		 push @new2,"*";
		}

	my @Replacement=@$ref_Replacement;
	my @Synonymous=@$ref_Synonymous;
	my @Damaging=@$ref_damaging;

	my $OrignalR= scalar(@Replacement);
	my $OrignalS= scalar(@Synonymous);
	#Silent first
	for (my $j=0;$j<scalar(@Synonymous);$j++)
		{
		 $Synonymous[$j]=$Synonymous[$j]-1;#make sure the position is exactly matching, need to subtract 1 for starting at 0
		 my $t= int($Synonymous[$j]/$scale);
		 $new2[$t]="S";
		}

#Damaging second
	for (my $j=0;$j<scalar(@Damaging);$j++)
		{
		 $Damaging[$j]=$Damaging[$j]-1;#make sure the position is exactly matching, need to subtract 1 for starting at 0
		 my $t= int($Damaging[$j]/$scale);
		 $new2[$t]="D";
		}
#Finally Replacement
	for (my $j=0;$j<scalar(@Replacement);$j++)
		{
		 $Replacement[$j]=$Replacement[$j]-1; #make sure the position is exactly matching, need to subtract 1 for starting at 0
		 my $t= int($Replacement[$j]/$scale);
		 $new2[$t]="R";
		}



	my $countR=0;
	my $countS=0;
	my $countD=0;


	for (my $ii=0;$ii<$length;$ii++)
		{
		 if ($new2[$ii] eq "R")
			 {
			  $countR=$countR+1;
			 }

		 if ($new2[$ii] eq "S")
			 {
			  $countS=$countS+1;
			 }

		 if ($new2[$ii] eq "D")
			 {
			  $countD=$countD+1;
			 }			 
		}

	printf "The total number of Replacement, Synonymous, Damaging after scaling is %d and %d.\n\n",$countR, $countS, $countD;
	printf ">%s_%s.%d_R_%d_S_%d_D_%d_len_%d_R_%d_S_%d\n%s\n",$geneID,$isoform,$length,$countR,$countS,$countD,$length1,$OrignalR,$OrignalS,join("",@new2);
	return @new2;
 
}
###############################################################################


###############################################################################
# make sure the position is exactly matching, need to subtract 1 for starting at 0
# $Replacement[$j]=$Replacement[$j]-1; 
sub CreateRS_ScaledSeq
{
	my ($scale,$geneID,$isoform,$length1,$ref_Synonymous,$ref_Replacement)=@_;

	my @new2;

	my $length=int($length1/$scale);
	for (my $i=0;$i<$length;$i++)
		{
		 push @new2,"*";
		}

	my @Replacement=@$ref_Replacement;
	my @Synonymous=@$ref_Synonymous;

	my $OrignalR= scalar(@Replacement);
	my $OrignalS= scalar(@Synonymous);
	for (my $j=0;$j<scalar(@Replacement);$j++)
		{
		 $Replacement[$j]=$Replacement[$j]-1; #make sure the position is exactly matching, need to subtract 1 for starting at 0
		 my $t= int($Replacement[$j]/$scale);
		 $new2[$t]="R";
		}

	for (my $j=0;$j<scalar(@Synonymous);$j++)
		{
		 $Synonymous[$j]=$Synonymous[$j]-1;#make sure the position is exactly matching, need to subtract 1 for starting at 0
		 my $t= int($Synonymous[$j]/$scale);
		 $new2[$t]="S";
		}

	my $countR=0;
	my $countS=0;

	for (my $ii=0;$ii<$length;$ii++)
		{
		 if ($new2[$ii] eq "R")
			 {
			  $countR=$countR+1;
			 }

		 if ($new2[$ii] eq "S")
			 {
			  $countS=$countS+1;
			 }
		}

	printf "The total number of Replacement and Synonymous after scaling is %d and %d.\n\n",$countR, $countS;
	printf ">%s_%s.%d_R_%d_S_%d_len_%d_R_%d_S_%d\n%s\n",$geneID,$isoform,$length,$countR,$countS,$length1,$OrignalR,$OrignalS,join("",@new2);
	return @new2;
 
}
###############################################################################


###############################################################################
sub PrintSeq{
my ($outdir,$ref_new,$GeneID,$length,$ID)=@_;
my @new_seq=@$ref_new;
my $nnew=join("",@new_seq);
my $count_replacement= ($nnew=~ tr/"R"//);
my $count_synonymous= ($nnew=~ tr/"S"//);
my $count_both= ($nnew=~ tr/"B"//);

my $total=$count_replacement+$count_synonymous+$count_both;
my $output_file_name=$GeneID."_".$length."_".$ID."_".$total.".txt";
my $contents=">".$GeneID.$length.$ID."total_".$total."_R_".$count_replacement."_S_".$count_synonymous."_B_".$count_both."\n".$nnew."\n";

    MAF::WriteFile($outdir,$output_file_name,$contents);
    
    return ($output_file_name);

}
###############################################################################





