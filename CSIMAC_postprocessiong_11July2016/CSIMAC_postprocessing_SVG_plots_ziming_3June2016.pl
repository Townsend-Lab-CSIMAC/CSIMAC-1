#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  SVG_plot_cMACPRF_ziming.pl
#
#        USAGE:  perl SVG_plot_cMACPRF_ziming.pl Combined_LUSC;
#                perl SVG_plot_cMACPRF_ziming.pl Combined_LUAD;
#                perl SVG_plot_cMACPRF_ziming.pl Gilead_LUSC;
#                perl SVG_plot_cMACPRF_ziming.pl Gilead_LUAD;
#        INPUT:  cMACPRF output
#       OUTPUT:  Plot with gamma, Confidence Intervals, and Mutation Status, silent -1, missense 1, recurrent 2.
#    ATTENTION:  Check the directory, $dir, $indir, $outdir; Check $TumorType and @GeneID.
#  DESCRIPTION:  Postprocess cMACPRF output to automate plot gamma values with 95% CI. The significant part will be marked as red or blue. The MA results will change from red to blue.
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  Fixed 11/05/2013 The pbs file, the recurrent file name - removing the directory. 
#        NOTES:  ---
#       AUTHOR:   (Zi-Ming Zhao), <ziming.gt@gmail.com>
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  April 24, 2013
#      UPDATED:  June 3rd, 2016; May 6th, 2016; November 13th, 2013; January 23, 2014
#===============================================================================

#REVISION RECORDS:
###################################################
#ZMZ 05/06/2016 revised sub-functions SVG_plot, Draw_info, Draw_x_recur, Draw_x_axis, and Draw_y_axis

#ZMZ 05/06/2016 Draw_y_axis, Label y axis, make it vertical by using transform=> 'translate(20,170) rotate(-90)'
#ZMZ 05/05/2016 SVG_plot, Re-Set the boundary for the plot
#ZMZ 05/05/2016 Draw_x_axis, move $Xlabel closer to x-axis, y=>($Yaxis_y2+25), 25 instead of 35
#ZMZ 05/05/2016 Draw_info: Revise plot for recurrent site gamma, and plot them separately from regional gamma, line from lower CI to upper CI, and blue dot the recurrent gamma.
#ZMZ 05/05/2016 Draw_info: Reassign start position for regional gamma plot, empty sites with recurrent gammas
#ZMZ 05/05/2016, in Draw_x_recur, y-5, not y+10, same as replacement mutations labels
#ZMZ 05/04/2016, in Draw_x_recur, position should just be $pos, not $pos+1 
#ZMZ 05/05/2016 Draw_x_axis, move $Xlabel closer to x-axis, y=>($Yaxis_y2+30), 30 instead of 35
#ZMZ 05/05/2016 Draw_y_axis, $Yvalues closer to y-axis, $Xaxis_x1-20, use 20 instead of 30.
###################################################

###################################################
#01/23/2014: change back not smooth the curve by skipping Non-value sites
# 01/17/2014 Draw_info, smooth the curve by skip Non-value sites, and use the neighboring sites with values to plot. Switch the site and data column in the for loop, to work on each data (gamma, 95% CI gamma) first for all sites, since Non-value sites need to be skipped. 
#? 11/08/2013 To improve: generate the plot for each gene separately, not to gather all the data for all genes first, then plot all at the same time.
#? 11/18/2013 To improve: get the hash with labels for the cMACPRF data, instead of using array of certain columns, since it is easy to be messed up with modified cMACPRF outputs.
# 11/07/2013 Updates: Revise to read the most updated cMACPRF outputs to extract the values of useful data; Draw recurrent mutations in the red color under the x axis.
# 11/08/2013 Fixed the y axis ticks and labels in the subroutine Draw_y_axis, x axis at the postion y=0, and have negative y axis ticks and labels.
# 11/09/2013 Fixed bug in the subroutine Draw_y_axis, while loop go up to the TickNum, since y=0 is always included, the criteria is also held by the boundary; Remove the drawing of y=0 outside the loop, since it is redundant.
# 11/11/2013 Minimum y --- use lower gamma CI, instead of model averaged gamma; use upper gamma CI for maximum y
# 11/13/2013 ?To Solve: don't plot with Null...Checked SVG plot Problems with Null: TMSB10 - Null, CAMK2N2 - Null, DCUN1D3 - Null, MT1H - Null;

# 11/13/2013 Todo: cMAC-PRF analyses of CBL and other top genes (as we have done for melanoma) for NSCLC?
# 11/15/2013 Fixed bug: y coordinates ticks, for tick between 1 and 10, if the digit is over 0.5, ground to the integer plus one. [gene TSKS as an example]

###################################################

###################################################
#April 24, 2013
#Max can be adjusted, min=0
#Generate several files.
#The scales for y axis should be integrate and for x axis better each step 50 or 100 ect.. 
#Add the nonsynonymous sites as dash to the x axis.
#Add the recurrent mutation numbers to the x axis. (recurrent mustaion is not using)
###################################################
#Note change directories from line 60 to line 92, and then line 286.

#!/usr/bin/perl -w
#use strict;
use SVG;
use Carp;
#use SVG::TT::Graph::BarLine;

##The directory information: input file directory, output file directory
#my $dir="/Users/zimingzhao/Desktop/cMAC-PRF/";


#GBM
#my $dir="/Users/gadarethhiggs/Desktop/Yale/Research/Rotations/R5-Fall 2014/CSIMAC/GBM_CSIMAC/";


#my $indir=$dir."GileadMelanoma/Gilead_Melanoma_outputs_01292014/";
#Melanoma
my $indir=$dir;
#GBM
#my $indir=$dir."GBM_CSIMAC_Results/Completed Genes/";
#my $outdirSVG=$dir."GileadMelanoma/cMACPRF_SVG_plots_Melanoma_01292014/";


#my $dir="/Users/zimingzhao/Dropbox/projects/csi-mac/CSIMAC_postprocessing_05022016/";
#my $outdir="/Users/zimingzhao/Dropbox/projects/csi-mac/CSIMAC_postprocessing_05022016/";
my $dir="/Users/zimingzhao/Dropbox/projects/CSI-MAC/CSIMAC_postprocessing_05022016/";
my $outdir="/Users/zimingzhao/Dropbox/projects/CSI-MAC/CSIMAC_postprocessing_05022016/";

#Melanoma
my @GeneIDs1;
my @GeneIDs2;
my @GeneIDs;

#GBM
#my $outdir=$dir."GBM_CSIMAC_SVG_plots/Top20/";

#Melanoma

#GBM (Glioblastoma)
#my $tumor_type="GBM";



#open(my $fh, "<", "CSIMAC_Gilead_Melanoma-top25.txt")
#MDAT Mutated Genes

#open(my $fh, "<", "MDAT_mutated_genes.txt")
#    or die "Failed to open file with gene IDs: $!\n";
#while(<$fh>) {
#    chomp;
#    push @GeneIDs1, $_;
#}
#close $fh;

#open(my $fh, "<", "complete_list.txt")
#    or die "Failed to open file with gene IDs: $!\n";
#while(<$fh>) {
#    chomp;
#    push @GeneIDs2, $_;
#}
#close $fh;

#my $tumor_type="lusc_pc";
my $tumor_type="brca_pc";
#my @GeneIDs = qw(TP53 AKT1 PIK3CA SF3B1);
#my @GeneIDs = qw(EGFR_regional);
my @GeneIDs = qw(PIK3CA);
#my @GeneIDs;
#open (my $fh,"<","complete_list.txt")
#or die "Failed to open file with gene IDs: $!\n";
#while(<$fh>) {
#    chomp;
#    push @GeneIDs, $_;
#}
#close $fh;



#print @GeneIDs;


my $scale=1;

#my $USAGE1 = "\nPlease, provide tumor type as inputs\n";
##my $tumor_type= $ARGV[0] or croak $USAGE1; # as part of the output
##my $indir=$dir.$tumor_type."_outputs_11212013/";

SVG_plot($tumor_type,\@GeneIDs);

######################################################################################### 
#ZMZ 05/05/2016 SVG_plot, Re-Set the boundary for the plot

sub SVG_plot
{
my ($TumorType,$ref_Genes)=@_;

my @GeneID=@$ref_Genes;

##Go over each gene, extract the cMACPRF data table, to save in the array_Files, and MutationStatus in array_Mut.
my @array_Files; ##The list of the cMACPRF data table files for each gene in the gene list
my @array_Mut; ##The list of the mutation status array reference for each gene in the gene list
my @array_Xticks;
my @array_max;
foreach my $ii (@GeneID) 
{ 
 my ($bin_size,$max_MA_r,$cMACPRF_data,@MutStatus)=CreateBICtable($ii,$TumorType,$indir);
 push @array_Files, $cMACPRF_data; ##Keep each gene's cMACPRF data file
 push @array_Mut, \@MutStatus;  ##Keep each gene's mutation status list's reference
 push @array_Xticks,$bin_size; ##Keep the bin size for each gene's cMACPRF plot
 push @array_max,$max_MA_r; ##Keep the maximum model averaged gamma for each gene
 }

########################################
##Go over each gene cMACPRF data file, to plot Model_Averaged Gamma, Lower 95% CI, Upper 95% CI, and the Mutation Status
for(my $i=0;$i<@array_Files;$i++){


    ###########################
    printf "#######################STEP1: Get the cMACPRF data file using hash_info\n";
    ##printf "The data file: %s\n",$array_Files[$i];
    ##Key:position of the gene
    ##Values:gamma/lower/upper/MutationStatus
    my %hash_info=();
    my $infile=$indir.$array_Files[$i];
    get_info($infile,\%hash_info);
    my $positions= scalar keys %hash_info;
    printf "##The number of sites in Gene %s: %d\n",$GeneID[$i],$positions;

    ###########################
    printf "#######################STEP2: Get mutation status: Replacement =1, Recurrent=2, Silent=-1.\n";
    ##Key: position (start from 0)
    ##Value:number
    my %hash_nonsyn=();
    my @GeneMutStatus;
    for my $each (@{$array_Mut[$i]})
    {
     push @GeneMutStatus,$each;
     #print "***",$each,"***\n";
    }

    my @Pos_nonsyn;
    my @Pos_recur;
    my $RepCount=0;
    my $RCurCount=0;
    
    #######################
    printf "#######################STEP3: Get the positions with mutation status\n";
	for(my $i=0;$i<scalar(@GeneMutStatus);$i++){
	    my $word=$GeneMutStatus[$i];
	    if($word eq '1'){
		push @Pos_nonsyn,$i;
		$RepCount++;
		#printf "Value: %s\t Count: %d\n",$word,$count;
	    }
		##Get the recurrent mutation information
		##Key: position (start from 0)                                                                                                                 
		##Value:number	    
	    if($word eq '2'){
		push @Pos_recur,$i;
		$RCurCount++;
		#printf "Value: %s\t Count: %d\n",$word,$count;
	    }
    }

    #######################
    printf "#######################STEP4: Creat a new SVG handle\n";
    my $svg=SVG->new(width =>600, height => 300);


    #######################
    ##Define the groups
    my $grp_axis= $svg->group(id => 'group_axis', style => {stroke=>'black',},);
    my $grp_data=$svg->group(id => 'group_data', style => {stroke=>'black'});
    my $grp_cis=$svg->group(id => 'group_cis', style => {stroke=>"rgb(0,255,127)", 'stroke-width'=>"0.5"},);
    my $grp_zero=$svg->group(id => 'group_zero', style => {stroke=>"rgb(0,255,0)",},);
    my $grp_text=$svg->group(id => 'group_text', style => {stroke=>'black'},);
    #Leave how much on the right side of x axis.
    my $x_gap=10;
    
    ##green color 0-255-127

    #######################
    printf "#######################STEP5: Draw the x & y axis lines, ticks, labels.\n";
    #ZMZ 05/05/2016 Set the boundary for the plot
    my $border_left=50; 
    my $border_right=590;
    my $border_up=10;
    my $border_down=260;
    my $Yaxis_x1=$border_left;
    my $Yaxis_y1=$border_up;
    my $Yaxis_x2=$border_left;
    my $Yaxis_y2=$border_down;
    my $Xaxis_x1=$border_left;
    my $Xaxis_y1=$border_down;
    my $Xaxis_x2=$border_right;
    my $Xaxis_y2=$border_down;
    my $tmp_x=$Xaxis_x2;
    my $tmp_y=$Yaxis_y1;
    ##Get the box area of the x and y
    ##my $tmp_x_axis=$grp_axis->line(id=>'l3',x1=>$Yaxis_x1,y1=>$Yaxis_y1,x2=>$tmp_x,y2=>$tmp_y); ##Up boundary of x axis
    ##my $tmp_y_axis=$grp_axis->line(id=>'l4',x1=>$tmp_x,y1=>$tmp_y,x2=>$Xaxis_x2,y2=>$Xaxis_y2); ##Right boundary of y axis

    #######################
    ###printf "#######################STEP: Get the X max and Y max & min\n";
    ##Get the X max which starts at 0.
    my @array=keys %hash_info;
    my $Xmax=@array-1;
    printf "##The maximum size of x axis sites: %d\n",$Xmax;

    ##Get the Y max & min
    
    my ($Ymax_r,$Ymin_r)=get_max_min_info(\%hash_info,0);
    my ($Ymax_l,$Ymin_l)=get_max_min_info(\%hash_info,1);
    my ($Ymax_u,$Ymin_u)=get_max_min_info(\%hash_info,2);
    printf "##Get gamma maximum %.2f, minimum %.2f\n",$Ymax_r, $Ymin_r;
    printf "##Get gamma lower CI maximum %.2f, minimum %.2f\n",$Ymax_l, $Ymin_l;
    printf "##Get gamma Upper CI maximum %.2f, minimum %.2f\n",$Ymax_u, $Ymin_u;

    my $Ymax=$Ymax_u; ##The maximum y axis is the maximum of MA gamma
    my $Ymin=$Ymin_l; ##The minimum y axis is the minimum of MA gamma
    my $Ydiff=$Ymax-$Ymin; ##The difference between the maximum and minimum gamma

    #######################
    ##printf "#######################STEP: Build the x & y axes\n";
    
    ##Move the x axis to the y=0
    my $Y_position0=Draw_y_axis($Ymax,$Ymin,$Yaxis_y1,$Yaxis_y2,$Xaxis_x1,$Xaxis_x2,$grp_axis,$svg,$positions);
    Draw_x_axis($Y_position0,$GeneID[$i],$Yaxis_y1,$Xmax,$Xaxis_x1,$Xaxis_x2,$Yaxis_y2,$x_gap,$grp_axis,$svg,$array_Xticks[$i]);
    
    my $yaxis=$grp_axis->line(id=>'l1',x1=>$Yaxis_x1,y1=>$Yaxis_y1,x2=>$Yaxis_x2,y2=>$Yaxis_y2); ##Left y axis from up to down
    my $xaxis=$grp_axis->line(id=>'l2',x1=>$Xaxis_x1,y1=>$Y_position0,x2=>$Xaxis_x2,y2=>$Y_position0); ##Down x axis from left to right

    #######################
    printf "#######################STEP6: Draw the mutation status onto x axis, Replacement and Recurrent\n";
    Draw_x_replacement($Xmax,$Xaxis_x1,$Xaxis_x2,$Y_position0,$x_gap,$svg,@Pos_nonsyn);
    Draw_x_recur($Xmax,$Xaxis_x1,$Xaxis_x2,$Y_position0,$x_gap,$svg,@Pos_recur);


    ########################
    printf "#######################STEP7: Draw the gamma value, and 95% CIs\n";
    Draw_info(\%hash_info,$Ymax_r,$Ymin_r,$Xmax,$Ymax,$Ymin,$Xaxis_x1,$Xaxis_x2,$Yaxis_y1,$Yaxis_y2,$x_gap,$Ydiff,$grp_cis,$grp_zero,$svg);

    ########################
    printf "#######################STEP8: Gene %s: Output the result into a svg file!\n",$GeneID[$i];
    my $out;
    my ($name)=($array_Files[$i]=~/(.*)_Table\.txt/);
    my $outfile=$outdir.$name.'.svg';
    open($out,">$outfile") or die "Can not open the file!\n";
    print $out $svg->xmlify;
    close($out);
}
exit;

}
######################################################################################### 




######################################################################################### 
##Extract the cMACPRF data from the cMACPRF output directly
sub CreateBICtable
{
my ($GeneID,$TumorType,$indir)=@_;
#Melanoma
#my $f=$indir.$GeneID.".".$TumorType."_bulldogn.CSIMAC.output.txt";
my $f=$indir.$GeneID.".".$TumorType.".CSIMAC.output.txt";
#GBM
    #my $f=$indir.$GeneID."_".$TumorType."_CSIMAC.output.txt";
printf "********Subroutine: Create BIC table for GeneID: $GeneID*********\n";
##printf "The cMACPRF output file as the input: %s\n",$f;

#######################
##Open the cMACPRF output file
if (!open(IN,$f)){
   croak ("Could not open file $f .\n")
}
#croak ('Could not open file $f .\n') unless open(IN,$f); 
my @c= <IN>;
close(IN);
my $c= join("",@c);

#######################
##Extract the data table from cMACPRF output file, keep only the columns of Positions, Gamma, LowerCI, UpperCI, mutationStatus
my $s= "Position\tGamma_Cancer\tLower_CI_Gamma_c\tUpper_CI_Gamma_c\tMutationStatus\n";

 my $loc1= index($c,"Position	      MS_DivRep	     MA_DivRep	     MS_DivSys	     MA_DivSys	     Gamma_Cancer   	Lower_CI_Gamma_c	Upper_CI_Gamma_c	MutationSymbol_Q2R1S-1*0	MutationStatus ");

 my $loc2= index($c,"\n",$loc1+2);
 my $loc3= index($c,"Abbreviation:", $loc1+2);
 my $newc="";
 if ($loc1!=-1 and $loc2!=-1 and $loc3!=-1) {
  $newc = substr($c,$loc1,$loc3-$loc1-1);
  }
  else { croak "Could not find the cMACPRF output in the right format in CreateBICtable! \n";}

$newc=~s/\n\n/\n/g; ##Remove the empty lines
$newc=~s/\n\n/\n/g; ##Remove the empty lines

my @c2=split(/\n/,$newc);
shift @c2; ##Remove the first line, the header line.
##printf "First item: ***%s***\n",$c2[0];
##printf "Last item: ***%s***\n",$c2[-1];

my $sep2="\t";
my $aa_length= scalar(@c2);
my $bin_size= int($aa_length/100)+1;
$bin_size=$bin_size*10*$scale; ##x axis bin size is the integer of amino acid length/10.

my @MA;
my @MutStatus;

##Go over each line of the cMACPRF output data to extract only the columns of interests, including position, MA gamma, CIs, mutationStatus.
for (my $j=0;$j<scalar(@c2);$j++)
{
 my @tt= split(/$sep2/,$c2[$j]);
 for my $each (@tt) { $each=~s/ //g;}
 my $temp= $tt[5]; ##Model averaged gamma value
 push @MA, $temp; ##Keep all model averaged gamma
 $tt[-2]=~s/ //g; ##Remove space
 push @MutStatus,$tt[8]; ##The mutation status column
 my $t1=$tt[0].$sep2.$tt[5].$sep2.$tt[6].$sep2.$tt[7].$sep2.$tt[8]; # Keep only the position[0], gamma[5], lower[6], upperCI[7], MutationStatus[8].
 $s.=$t1."\n";
 ###printf "*Pos: *%s*; *Gamma: *%s*; *LCI: *%s*; *UCI: *%s*; *Status: *%s******\n",$tt[0],$tt[5],$tt[6],$tt[7],$tt[8];
 
}

#######################
##Sort gamma values, and keep the maximum gamma value
my $max_MA_r= 0;
my @new_MA= sort {$b <=> $a} @MA;
$max_MA_r= $new_MA[0];
printf "##The gene size, bin size, maximum and minimum MA gamma: %d, %d, %.2f,%.2f.\n",$aa_length,$bin_size,$new_MA[0],$new_MA[-1];
#printf "Sorted MA: ", (join "\t",@new_MA), "\n";
if ($max_MA_r<0.05) {$max_MA_r=0.05;}

#######################
##Save the extracted interested gamma data in the new file
my $new_f=$GeneID.".".$TumorType.'_cMACPRF_Table.txt';
my $new=">".$indir.$new_f;
open (OUT, $new);
my $pres= select(OUT);
print $s;
select $pres;
close(OUT);

##Return the following items
return ($bin_size,$max_MA_r,$new_f,@MutStatus);
 
}
######################################################################################### 


######################################################################################### 
sub get_info{
    my ($file,$H)=@_;
    
    #######################
    open(FH,"$file") or die "Can not open the file $file:$!\n";
    ##printf "********Subroutine: Get information from file: %s\n",$file;
    my $count=0;
    my $title='';
    
    #######################
    while(my $line=<FH>){
        $count++;
        
        chomp($line);
        my @array=split(/\t/,$line);
        #==========================
        for(my $i=0;$i<@array;$i++){
            $array[$i]=~s/\s|\t//g;
        }
        #==========================

        if($count==1){
            $title=$array[-4];
            next;
        }
        
        #######################
        ##Keep in the hash the desired values
        ##Gamma/lower/upper/MutationStatus
        $$H{$count-2}[0]=$array[1];
        $$H{$count-2}[1]=$array[2];
        $$H{$count-2}[2]=$array[3];
        $$H{$count-2}[3]=$array[4];
        
##GAHEdit Error checking
#       printf "The output of the first array is: $$H{count-2}[0]";
##GAHEdit Error checking
        
    }
    close(FH);
    return;
}
######################################################################################### 


######################################################################################### 
##Get the maximum and minimum values for the data, for getting the proper tick of the values
##my ($Ymax,$Ymin)=get_max_min_info(\%hash_info,0);
sub get_max_min_info{
    my ($H,$num)=@_;

    #######################
    ##Find the Max
    my $Ymax = -100000;
    foreach my $key (keys %{$H}){
	if($$H{$key}[$num] ne 'NULL' and $$H{$key}[$num] ne 'N-INF' and $$H{$key}[$num] ne 'INF'){
	    if($$H{$key}[$num]>$Ymax){
		$Ymax=$$H{$key}[$num];
	    }
	}
    }

    #######################
    ##Find the min
    my $Ymin=100000;
    foreach my $key (keys %{$H}){
	if($$H{$key}[$num] ne 'NULL' and $$H{$key}[$num] ne 'N-INF' and $$H{$key}[$num] ne 'INF'){
	    if($$H{$key}[$num]<$Ymin){
		$Ymin=$$H{$key}[$num];
	    }
	}
    }
    ##printf "********Subroutine: Get the maximum and minimum value: %.2f,%.2f\n",$Ymax,$Ymin;
    
    return($Ymax,$Ymin);
}
######################################################################################### 


#########################################################################################    
##Draw x axis, the position, the ticks, the maximum x, the label of x axis 
#ZMZ 05/05/2016 Draw_x_axis, move $Xlabel closer to x-axis, y=>($Yaxis_y2+25), 25 instead of 35
sub Draw_x_axis{
    printf "********Subroutine: Draw x axis\n";
    my ($YPosition0,$gene,$Yaxis_y1,$Xmax,$Xaxis_x1,$Xaxis_x2,$Yaxis_y2,$x_gap,$grp_axis,$svg,$Xticks)=@_;

    #######################
    ##Get the value of every scale of the X axis
    $Xmax=$Xmax*$scale;
    my $times=int($Xmax/$Xticks);

    #######################
    ##Draw the scale of the X axis       
    for(my $i = 0 ; $i <= $times; $i++)
    {
	my $Xposition = $Xaxis_x1 + ($i * $Xticks * ($Xaxis_x2-$Xaxis_x1-$x_gap)) / ($Xmax - 1);
	my $tmp_tick=$grp_axis->line(x1=>$Xposition,y1=>($YPosition0+5),x2=>($Xposition),y2=>($YPosition0),);
	#printf "Tmp_tick in Draw x axis: %s\n",$tmp_tick;
	my $text=0;
	if($i==0){
	    $text=$svg->text(x=>$Xposition+2,y=>($YPosition0+15),style => {'font-size'=>'10'})->cdata($i*$Xticks);
        }else{
            $text=$svg->text(x=>($Xposition-9), y=>($YPosition0+15),style => {'font-size'=>'10'})->cdata($i*$Xticks);
        }
     }
     
     #######################
     my $Xlabel=$gene." Coordinates (aa)";
     my $Xposition1 =$Xaxis_x1 + ($times/5*2 * $Xticks * ($Xaxis_x2-$Xaxis_x1-$x_gap)) / ($Xmax - 1);
     my $text1=$svg->text(x=>$Xposition1,y=>($Yaxis_y2+30),style => {'font-size'=>'10'})->cdata($Xlabel);

     return;
}
#########################################################################################

######################################################################################### 
#??? How to label the y=0, and both positive and negative based on the y values
##??? Draw the x axis on the y=0
##Draw y axis, the position, the ticks, the maximum y, the label of y axis
#ZMZ 05/05/2016 Draw_y_axis, $Yvalues closer to y-axis, $Xaxis_x1-20, use 20 instead of 30.
###ZMZ 05/05/2016 Draw_y_axis, To Label y axis

sub Draw_y_axis{
    my ($Ymax,$Ymin,$Yaxis_y1,$Yaxis_y2,$Xaxis_x1,$Xaxis_x2,$grp_axis,$svg,$sites)=@_;
    printf "********Subroutine: Draw y axis ticks and labels.\n";
    
    #######################
    ##printf "Get the value of the every scale of the Y axis\n";
    my $TickNum=10; ##the number of ticks on the y axis

    $Ymax=sprintf "%.2f",$Ymax; ##Keep only two digits of maximum Y
    $Ymin=sprintf "%.2f",$Ymin; ##Keep only two digits of minimum Y
    ##printf "Minimum and Maximum Y: %.2f, %.2f\n",$Ymin,$Ymax;
    my $Ydiff=$Ymax;
    my $Yposition0=$Yaxis_y2; ##If y min is greater than 0, then y starts from coordinates 0, ydiff is Ymax.

    ##If y min is less than 0, then the y=0 is between Yaxis_y2 and Yaxis_y1
    if ($Ymin<0) {
        $Ydiff=$Ymax-$Ymin;
        $Yposition0 = $Yaxis_y2 + int(($Ymin* ($Yaxis_y2-$Yaxis_y1)) / $Ydiff);
        ##printf "Ymin less than 0! The position for y=0: %.2f\n",$Yposition0; 
    }
    if($Ydiff==0){
        $Ydiff=0.0001;
    }
    my $Yticks=$Ydiff/$TickNum; ##The length of each bar
    if ($Yticks>=10) {$Yticks=int($Yticks);}
    ##To make sure Ytick is 3 if Ytick=2.6, since if it grounds to 2, then the y axis label is not complete; get the big integer if it pass .5
    if ($Yticks>1 and $Yticks<10) { 
        my $Yticks_tmp=int($Yticks); 
        if ($Yticks-$Yticks_tmp>=0.5){
            $Yticks=$Yticks_tmp+1;
        }
        else {
            $Yticks=int($Yticks);
        }

    }


    #######################
    ##printf "Draw the scale of the Y axis\n";
    my $k=0;
    my $y_var=$Yaxis_y2-$Yaxis_y1+1; ##The length in the plot on the left y axis
    my $Yvalues=0;
    ##my $x_tick=($Xaxis_x2-$Xaxis_x1)/$sites; 
    ##printf "The length of one site in the plot on the x axis: %.2f\n",$x_tick;

    ##Draw the ticks on the y axis from y=0 to top till hitting the upper boundary of the plot; for negative Ymin, additional drawing from the y=0 position to the bottom till hitting the lower boundary of the plot;
    ##Add the corresponding number on the y axis ticks
    
    ###ZMZ 05/06/2016 Draw_y_axis, Label y axis, make it vertical by using transform=> 'translate(20,170) rotate(-90)'
    my $Ylabel="Selection Intensity";
    $text = $svg->text(transform=> 'translate(20,170) rotate(-90)', style => {'font-size'=>'10'})->cdata($Ylabel);

    

    ##Draw ticks for all Ymin>=0 and the ticks above 0 for Ymin<0; no ticks above the upper boundary of Yaxis_y1
    while($y_var>=$Yaxis_y1 and $k<=$TickNum){
        my $Yposition = $Yposition0 - int(($k * $Yticks * ($Yaxis_y2-$Yaxis_y1)) / $Ydiff);
        $y_var=$Yposition; ##The position for the current tick
        ##Quit and get out of the while loop if the tick is out of the upper boundary
        if($y_var<$Yaxis_y1){
            last;
        }

        my $tmp_tick=0;
        $tmp_tick=$grp_axis->line(x1=>$Xaxis_x1,y1=>$Yposition,x2=>($Xaxis_x1-5),y2=>$Yposition); ##Draw the ticks with the width of 5 on the y axis
	    $Yvalues=$Yticks*$k;##The value for each tick above 0
        $Yvalues=sprintf "%.2f",$Yvalues; ##Keep only two digits of numbers for the y axis labels
        if ($Yticks>1) {$Yvalues=sprintf "%d",$Yvalues; }
        my $text='';
        $text=$svg->text(x=>($Xaxis_x1-20),y=>($Yposition+2),style => {'font-size'=>'10'})->cdata($Yvalues); ##Label the text next to the y axis ticks
        $k++; ##Move to the next tick above the current positive one
        printf "The ticks: %d, y axis: %.2f\n",$k,$Yvalues;
     }
    
    ##Draw ticks for Ymin<0, the ticks below y=0; no  ticks below the lower bound Yaxis_y2
    while($Ymin<0 and $y_var<=$Yaxis_y2 and $k<=$TickNum){
        my $Yposition = $Yposition0 + int(($k * $Yticks * ($Yaxis_y2-$Yaxis_y1)) / $Ydiff);
        #If this y bar is out of the figure, we stop drawing.
        $y_var=$Yposition;
        ##Quit and get out of the while loop if the ticks move out of the lower bound
        if($y_var>$Yaxis_y2){
            last;
        }

        my $tmp_tick=0;
        $tmp_tick=$grp_axis->line(x1=>$Xaxis_x1,y1=>$Yposition,x2=>($Xaxis_x1-5),y2=>$Yposition); ##Draw the tick's line with the width of 5 on the y axis
        $Yvalues=-$Yticks*$k; ##The value for each tick below 0
        $Yvalues=sprintf "%.2f",$Yvalues; ##Keep only two digits of numbers for the y axis labels
        if ($Yticks>1) {$Yvalues=sprintf "%d",$Yvalues; }
        my $text='';
        $text=$svg->text(x=>($Xaxis_x1-20),y=>($Yposition+2),style => {'font-size'=>'10'})->cdata($Yvalues); ##Label the text next to the y axis ticks
        $k++; ##Move to the next tick, below the 0 and current negative one
        ###printf "The ticks: %d, y axis: %.2f\n",$k,$Yvalues;
     }   
     #my $yTitle=$svg->secondary_y_title(x=>580,y=>($Yaxis_y2),style => {'font-size'=>'12'})->cdata('Selection Intensity by cMACPRF');
     #printf "$y_var\n";
     return ($Yposition0);
}
######################################################################################### 



######################################################################################### 
##Draw the positions of replacement mutations on the plot with the blue color rgb(0,0,255), above the x axis, with the height of 5
##Draw_x_replacement($Xaxis_x1,$Xaxis_x2,$Yaxis_y2,$svg,\%hash_nonsyn);
sub Draw_x_replacement{
printf "********Subroutine: Draw_x_replacement\n";
    my ($Xmax,$Xaxis_x1,$Xaxis_x2,$Yaxis_y2,$x_gap,$svg,@HH)=@_;

    foreach my $pos (@HH) {
   
	my $x=$Xaxis_x1 + (($pos+1)*($Xaxis_x2-$Xaxis_x1-$x_gap)) / ($Xmax - 1); ##Get the mutation site position on the x axis
	#printf "Mutaiton position: %d\t X position: %.2f\n",$pos+1,$x;
	my $tmp_line=$svg->line(x1=>$x,y1=>$Yaxis_y2,x2=>$x,y2=>($Yaxis_y2-5),style => {stroke=>"rgb(0,0,255)",'stroke-width'=>"1"},); # In blue color rgb(0,0,255)

     }                                                                                                                                                       
     return;                                                                                                                                                 
} 
#########################################################################################


#########################################################################################       
##Draw the positions of recurrent mutations on the plot with the red color rgb(255,0,0), below the x axis, with the height of 10
#ZMZ 05/04/2016, in Draw_x_recur, position should just be $pos, not $pos+1 
#ZMZ 05/05/2016, in Draw_x_recur, y-5, not y+10, same as replacement mutations labels
sub Draw_x_recur{
printf "********Subroutine: Draw_x_recurrent:\n";
    my ($Xmax,$Xaxis_x1,$Xaxis_x2,$Yaxis_y2,$x_gap,$svg,@HH)=@_;

    foreach my $pos (@HH) {
    #ZMZ 05/04/2016, position should just be $pos, not $pos+1
	my $x=$Xaxis_x1 + (($pos)*($Xaxis_x2-$Xaxis_x1-$x_gap)) / ($Xmax - 1); ##Get the mutation site position on the x axis
	#printf "Mutaiton position: %d\t X position: %.2f\n",$pos+1,$x;
	#ZMZ 05/05/2016, in Draw_x_recur, y-5, not y+10, same as replacement mutations labels
	my $tmp_line=$svg->line(x1=>$x,y1=>$Yaxis_y2,x2=>$x,y2=>($Yaxis_y2-5),style => {stroke=>"rgb(0,0,255)",'stroke-width'=>"1"},);

     }                                                                                                                                                       
     return;                                                                                                                                                 
}        
#########################################################################################

##Draw_info(\%hash_info,$Ymax_r,$Ymin_r,$Xmax,$Ymax,$Ymin,$Xaxis_x1,$Xaxis_x2,$Yaxis_y1,$Yaxis_y2,$x_gap,$Ydiff,$grp_cis,$grp_zero,$svg);
#ZMZ 05/05/2016 Draw_info: Revise plot for recurrent site gamma, and plot them separately from regional gamma, line from lower CI to upper CI, and blue dot the recurrent gamma.
#ZMZ 05/05/2016 Draw_info: Reassign start position for regional gamma plot, empty sites with recurrent gammas
######################################################################################### 
sub Draw_info{
my ($H,$Ymax_r,$Ymin_r,$Xmax,$Ymax,$Ymin,$Xaxis_x1,$Xaxis_x2,$Yaxis_y1,$Yaxis_y2,$x_gap,$Ydiff,$grp_cis,$grp_zero,$svg)=@_;
my $NonValueStarts=0;
my $NonValueEnds=0;
my $start=0;
my $end=0;
#######################
##Draw each data
#j is for each column of data, particularly gamma and 95% CI gamma;
#Go over each site in one dataset for gamma, then 95% CI gamma
#ZMZ 05/05/2016 create a list of y-location for lower CI gamma to plot starting from lower CI gamma
my @lowerCI_y;
my $k=0;
for(my $j=0;$j<=@{$$H{0}};$j++){
	printf "\nThe column of the data (gamma & 95% CIr): %d\n",$j;
	$k=0;
    #i is for each site
    $start=0;
    for(my $i=1;$i<=$Xmax;$i++){
        $start =$i-1;
        #ZMZ 05/04/2016 Draw_info: Revise plot for recurrent site gamma, and plot them separately from regional gamma.
        if ($$H{$i}[3]==2) {
			##printf "Recurrent Position %d: %d\n",$i,$$H{$i}[3];        
			my $x=$Xaxis_x1 + ($i * ($Xaxis_x2-$Xaxis_x1-$x_gap)) / ($Xmax - 1);
			my $y2=$Yaxis_y2 - (($$H{$i}[$j] - $Ymin) * ($Yaxis_y2-$Yaxis_y1)) / $Ydiff;
			#Draw 95% CI gamma - the green bar
			if ($j==1 or $j==2){my $tmp_line=$svg->line(x1=>$x-1,y1=>$y2,x2=>$x+1,y2=>$y2,style => {stroke=>"rgb(0,255,127)",'stroke-width'=>"1"},);}
			   
			#Keep the record of y axis position for the lower gamma
			if ($j==1) { push @lowerCI_y, $y2; } 
			#Draw the vertical line from lower CI gamma to higher CI gamma
			if ($j==2){ my $tmp_line=$svg->line(x1=>$x,y1=>$lowerCI_y[$k],x2=>$x,y2=>$y2,style => {stroke=>"rgb(0,255,127)",'stroke-width'=>"0.5"},); }
			$k++;
			next; 
        }

        #ZMZ 05/04/2016 Draw_info: Reassign start position for regional gamma plot, by checking up to the previous position for the status of recurrent
        if ($$H{$i}[3]!=2 and $$H{$i-1}[3]==2) {
        $start =$i;        
        }                                
        $end=$i;
        
		###GAHEdit Error Checking
			#      if($j==0){
			#        printf "$$H{$start}[$j]";
			#    }
		###GAHEdit Error Checking

		my $x1=$Xaxis_x1 + ($start * ($Xaxis_x2-$Xaxis_x1-$x_gap)) / ($Xmax - 1);
		my $x2=$Xaxis_x1 + ($end * ($Xaxis_x2-$Xaxis_x1-$x_gap)) / ($Xmax - 1);
		#printf "Position: %d; X axis: %.2f\n",$i-1,$x1;

	    my $y1=$Yaxis_y2 - (($$H{$start}[$j] - $Ymin) * ($Yaxis_y2-$Yaxis_y1)) / $Ydiff;
	    my $y2=$Yaxis_y2 - (($$H{$end}[$j] - $Ymin) * ($Yaxis_y2-$Yaxis_y1)) / $Ydiff;
        
        #Plot regional gamma 
        if($j==0){	    
			#######################
			#my $tmp_line=$svg->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,style => {stroke=>"rgb($r,$g,$b)",'stroke-width'=>"2"},);
			#ZMZ 05/05/2016 always make regional gamma the blue color
			if ($$H{$i}[3]!=2 and $$H{$i-1}[3]!=2) {
			my $tmp_line=$svg->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,style => {stroke=>"rgb(0,0,255)",'stroke-width'=>"1"},);
			#printf "Draw regional gamma!\n";
			}
			#ZMZ 05/05/2016 for regional gamma, empty the space if the previous one is recurrent, draw a flat bar for the current non-recurrent site
			if ($$H{$i}[3]!=2 and $$H{$i-1}[3]==2) { my $tmp_line=$svg->line(x1=>$x1-0.5,y1=>$y1,x2=>$x2+0.5,y2=>$y2,style => {stroke=>"rgb(0,0,255)",'stroke-width'=>"1"},); }  
        }
        #Plot 95% CI for regional gamma 
        elsif($j==1 or $j==2){
			if($$H{$start}[$j]>$Ymax){
				$y1=$Yaxis_y2 - (($Ymax - $Ymin) * ($Yaxis_y2-$Yaxis_y1)) / $Ydiff;
			}
			if($$H{$end}[$j]>$Ymax){
				$y2=$Yaxis_y2 - (($Ymax - $Ymin) * ($Yaxis_y2-$Yaxis_y1)) / $Ydiff;
			}
			my $tmp_line=$grp_cis->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2);		
	    }
	}
}
#ZMZ 05/05/2016 Draw_info: Draw recurrent site gamma at the last, to make it to the front
for(my $i=1;$i<=$Xmax;$i++){
	if ($$H{$i}[3]==2) {
		##printf "Recurrent Position %d: %d\n",$i,$$H{$i}[3];        
		my $x=$Xaxis_x1 + ($i * ($Xaxis_x2-$Xaxis_x1-$x_gap)) / ($Xmax - 1);
		my $y2=$Yaxis_y2 - (($$H{$i}[$j] - $Ymin) * ($Yaxis_y2-$Yaxis_y1)) / $Ydiff;
		if ($j==0){ my $tmp_line=$svg->circle(cx=>$x,cy=>$y2,r=>0.5, style => {stroke=>"rgb(0,0,255)"},);} 
	} 
}  
return;
}
######################################################################################### 

