#
#===============================================================================
#
#         FILE:  MAF_ziming.pm
#
#  DESCRIPTION:  This modules is used for parsing MAF (mutation annotation
#                files) files from TAGC database
#
#        FILES:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:   (Zi-Ming Zhao), <ziming.gt@gmail.com>
#      COMPANY:  
#      VERSION:  0.3
#      CREATED:  01/16/2013 
#      UPDATED:  02/17/2014
#     REVISION:  ---
#===============================================================================

# 02/17/2-14: For Replacement, considering other mutation types in addition to Missense, change codes accordingly in ExtractMutPos; 
# $interesting_muts_rx = qr/Silent|Missense|Nonsense|Frame_Shift|In_Frame/xmi;
# my $missense_rx = qr/Missense|Nonsense|Frame_Shift|In_Frame/xmi;
# my $virant_classification_rx = qr/SNP|INS|DEL/xmi; 
# 02/18/2014 $first_line =~ s/\r\n/\n/g; #The problem of Windows (CRLF) \r\n issue 

package MAF;

#use strict;
use warnings;
use Carp;
use Array::Utils;
use Data::Dumper;

=begin GHOSTCODE
=cut

############## CHECK "@default_labels33" FOR Hugo_Symbol ##############

#Gilead database TCGA MAF labels column labels, total 20
my @default_labels33 = qw(Hugo_Symbol	Entrez_Gene_Id	NCBI_Build	Chrom	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	dbSNP_RS	dbSNP_Val_Status	Tumor_Sample_Barcode	Normal_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2	Verification_Status	Validation_Status);
#Gildead database yale MAF labels column labels, total 18
my @default_labels3=qw(Hugo_Symbol	Entrez_Gene_Id	NCBI_Build	Chrom	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	dbSNP_RS	dbSNP_Val_Status	Tumor_Sample_Barcode	Normal_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2);

# ovary cancer label - panCancer labels
my @default_labels333=qw(Hugo_Symbol	Entrez_Gene_Id	Center	NCBI_Build	Chromosome	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	dbSNP_RS	dbSNP_Val_Status	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2	Tumor_Validation_Allele1	Tumor_Validation_Allele2	Match_Norm_Validation_Allele1	Match_Norm_Validation_Allele2	Verification_Status	Validation_Status	Mutation_Status	Sequencing_Phase	Sequence_Source	Validation_Method	Score	BAM_File	Sequencer	chromosome_name	start	stop	reference	variant	type	gene_name	transcript_name	transcript_species	transcript_source	transcript_version	strand	transcript_status	trv_type	c_position	amino_acid_change	ucsc_cons	domain	all_domains	deletion_substructures	transcript_error);

my @default_labels=@default_labels333;
#my @default_label_meningionmas=qw(Hugo_Symbol Entrez_Gene_Id  NCBI_Build  Chrom  Start_position  End_position  Strand  Variant_Classification  Variant_Type  Reference_Allele  Tumor_Seq_Allele1 Tumor_Seq_Allele2 dbSNP_RS  dbSNP_Val_Status  Tumor_Sample_Barcode  Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1  Match_Norm_Seq_Allele2);
#my @default_labels3=@default_label_meningionmas;
#RefFlat labels

# GBM MAX cancer label - panCancer labels (xinru:10.31.2015)
my @default_labels3333=qw(Hugo_Symbol	geneid	chr	start	strand	ref_allele	newbase	Variant_Classification	effect	patient);
# Hugo_Symbol	geneid	chr	start	strand	ref_allele	newbase	Variant_Classification	effect	patient


my @ref_labels=qw(geneName	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds);

my @MutsigcvLabels= qw(gene	N_nonsilent	N_silent	N_noncoding	n_nonsilent	n_silent	n_noncoding	nnei	x	X	r_x_X	r_n_N	r_nns_Nns	r_ns_Ns	r_nnc_Nnc	r_nsnc_Nsnc);
###############################################################################


#####################################
# Parsing a whole file - PanCancer database ??? columns
# and returns a ref of hash of hashes
#####################################

sub parse_MAF3333 {
    #Only input MAF filename
    my ($inmaf) = @_;
    
    # Specify field separator...
    my $RECORD_SEPARATOR = qq{\t};
    ### Need to change to the right number of columns 
    my $FIELD_COUNT      = 10;

    
    open my $MAF, '<', "$inmaf"
        or croak "$0: failed to open infile '$inmaf': $!\n";

    #checking labels of infile
    my $first_line = <$MAF>;
    $first_line =~ s/\r\n/\n/g; #The problem of Windows (CRLF) \r\n issue
    chomp $first_line;
    printf "First line: %s***\n", $first_line;
    my @labels = split $RECORD_SEPARATOR, $first_line, $FIELD_COUNT+1;
    my $size=scalar(@labels);
    my $i=0;
        for my $each (@labels)
    {
    $i++;
    $each =~ s/ //g;
    printf "**Label %d: %s, %s^^^\n",$i,$each,$default_labels3333[$i-1];
    }
        
    #printf "Last item:%s****\n",$labels[$size-1];
    #printf "Fourth item: ***%s****\n",$labels[3];
    #$labels[4] =~ s/Chromosome/Chrom/g;

    my $labels_num = scalar @labels;
    printf "Label number: %d; Field Count: %d\n",$labels_num,$FIELD_COUNT;
    croak "*******Different number of labels for [$inmaf]*****\n"
        if $labels_num != $FIELD_COUNT;
    
    #print "\nFIELDS: [@labels]\n";
    print "Total fields ::: $labels_num\n";

    my @new_labels = Array::Utils::array_minus(@labels, @default_labels3333);
    printf "The length of new labels, %d, %s\n",scalar(@new_labels), @new_labels;

    croak "@new_labels labels are not default labels" if @new_labels;
    
    
    my @mutations;
    my $record_num = 0;
    
    ## GAH Edit: Add list called geneIDs
    my @geneids;
    my $i=0;

    while (my $record = <$MAF>) {
        $i++;
        next if $record !~ /\S/;
        $record=~ s/\r\n/\n/g;
        chomp $record;
        #split columns
        my @values = split $RECORD_SEPARATOR, $record, $FIELD_COUNT+1;
        
        #removing leading and trailing white spaces
        s{^\s+|\s+$}{}g foreach @values;

        #Error if there is no Hugo symbol for that gene
        #croak "No hugo symbol for record num [$record_num]: $record"
        #    if !$values[0];
        #Error if there is no Hugo symbol for that gene

		if (!$values[0]) 
		{
		 print "No hugo symbol for record num [$record_num]: $record\n";
		 next;
		}
        else {
            push @geneids, $values[0];
        }
        #printf "line: %d, Gene name: %s\n",$i, $values[0];
        
        
        
        #renameing chr X and Y
        $values[2] =~ s/chromosome//i;
        $values[2] =~ s/chr//i;
        $values[2] = uc $values[2] ;
        $values[2] = 24 if $values[2] =~ /x/i;
        $values[2] = 25 if $values[2] =~ /y/i;
        $values[2] = 23 if $values[2] =~ /MT/i;
        $values[2] = 23 if $values[2] =~ /Mito/i;

        # chromosome index needs to be consistent with "UCSC_refseq_refFlat_hg19" file (chrX)
        # Mito not found in "UCSC_refseq_refFlat_hg19"
        $values[2] = 'X' if $values[2] =~ /24/;
        $values[2] = 'Y' if $values[2] =~ /25/;
        $values[2] = 'Mito' if $values[2] =~ /23/; 
        #printf "line: %d, chromosome %d\n",$i, $values[4];
        
        #some mutations are named missense_mutations and other just missense
        
        $values[7] =~ s/_Mutation$//xm;
        #printf "line: %d, Mutations: *%s***\n", $i,$values[8];
        #printf "Value: %s~~~~~~\n", join "\t", @values;
        
        #Tumor sample name - changing to Patient ID
        ###my %barcode_value_for = parse_barcode($values[14]);
        ###$values[14] = $barcode_value_for{PatientID};

        my $values_num = scalar @values;
        my %record_info;
        @record_info{@default_labels3333} = @values;
        #print "*******",join "*\t^", @values, "*****\n";

        #Error if we found a different number of columns
        croak "[$values_num] different from expected [$FIELD_COUNT]" 
            if $values_num != $FIELD_COUNT;
                
        #Populating a hash of array of hashes with records. 
        push @mutations, \%record_info;
        ++$record_num;

    } #end of while loop

    my $sizeGenes=scalar @mutations;
    printf "The total number of genes extracted from the file %s is %d.\n",$inmaf, $sizeGenes;
    print "\nFor file [$inmaf] we found $record_num mutations\n\n\n";

    close $MAF
        or carp "$0: failed to close infile '$inmaf': $!\n";

    return (\@mutations, \@geneids);
} 
###################################

########CHECK parse_MAF333 for

sub parse_MAF333 {
    #Only input MAF filename
    my ($inmaf) = @_;
    
    # Specify field separator...
    my $RECORD_SEPARATOR = qq{\t};
    ### Need to change to the right number of columns 
    my $FIELD_COUNT      = 53;

    
    open my $MAF, '<', "$inmaf"
        or croak "$0: failed to open infile '$inmaf': $!\n";

    #checking labels of infile
    my $first_line = <$MAF>;
    $first_line =~ s/\r\n/\n/g; #The problem of Windows (CRLF) \r\n issue
    chomp $first_line;
    printf "First line: %s***\n", $first_line;
    my @labels = split $RECORD_SEPARATOR, $first_line, $FIELD_COUNT+1;
    my $size=scalar(@labels);
    my $i=0;
        for my $each (@labels)
    {
    $i++;
    $each =~ s/ //g;
    printf "**Label %d: %s, %s^^^\n",$i,$each,$default_labels333[$i-1];
    }
        
    #printf "Last item:%s****\n",$labels[$size-1];
    #printf "Fourth item: ***%s****\n",$labels[3];
    #$labels[4] =~ s/Chromosome/Chrom/g;

    my $labels_num = scalar @labels;
    printf "Label number: %d; Field Count: %d\n",$labels_num,$FIELD_COUNT;
    croak "*******Different number of labels for [$inmaf]*****\n"
        if $labels_num != $FIELD_COUNT;
    
    #print "\nFIELDS: [@labels]\n";
    print "Total fields ::: $labels_num\n";

    my @new_labels = Array::Utils::array_minus(@labels, @default_labels333);
    printf "The length of new labels, %d, %s\n",scalar(@new_labels), @new_labels;

    croak "@new_labels labels are not default labels" if @new_labels;
    
    
    my @mutations;
    my $record_num = 0;
    
    ## GAH Edit: Add list called geneIDs
    my @geneids;
    my $i=0;

    while (my $record = <$MAF>) {
        $i++;
        next if $record !~ /\S/;
        $record=~ s/\r\n/\n/g;
        chomp $record;
        #split columns
        my @values = split $RECORD_SEPARATOR, $record, $FIELD_COUNT+1;
        
        #removing leading and trailing white spaces
        s{^\s+|\s+$}{}g foreach @values;

        #Error if there is no Hugo symbol for that gene
        #croak "No hugo symbol for record num [$record_num]: $record"
        #    if !$values[0];
        #Error if there is no Hugo symbol for that gene

		if (!$values[0]) 
		{
		 print "No hugo symbol for record num [$record_num]: $record\n";
		 next;
		}
        else {
            push @geneids, $values[0];
        }
        #printf "line: %d, Gene name: %s\n",$i, $values[0];
        
        
        
        #renameing chr X and Y
        $values[4] =~ s/chromosome//i;
        $values[4] =~ s/chr//i;
        $values[4] = uc $values[4] ;
        $values[4] = 24 if $values[4] =~ /x/i;
        $values[4] = 25 if $values[4] =~ /y/i;
        $values[4] = 23 if $values[4] =~ /MT/i;
        $values[4] = 23 if $values[4] =~ /Mito/i;

        # chromosome index needs to be consistent with "UCSC_refseq_refFlat_hg19" file (chrX)
        # Mito not found in "UCSC_refseq_refFlat_hg19"
        $values[4] = 'X' if $values[4] =~ /24/;
        $values[4] = 'Y' if $values[4] =~ /25/;
        $values[4] = 'Mito' if $values[4] =~ /23/; 
        #printf "line: %d, chromosome %d\n",$i, $values[4];
        
        #some mutations are named missense_mutations and other just missense
        
        $values[8] =~ s/_Mutation$//xm;
        #printf "line: %d, Mutations: *%s***\n", $i,$values[8];
        #printf "Value: %s~~~~~~\n", join "\t", @values;
        
        #Tumor sample name - changing to Patient ID
        ###my %barcode_value_for = parse_barcode($values[14]);
        ###$values[14] = $barcode_value_for{PatientID};

        my $values_num = scalar @values;
        my %record_info;
        @record_info{@default_labels333} = @values;
        #print "*******",join "*\t^", @values, "*****\n";

        #Error if we found a different number of columns
        croak "[$values_num] different from expected [$FIELD_COUNT]" 
            if $values_num != $FIELD_COUNT;
                
        #Populating a hash of array of hashes with records. 
        push @mutations, \%record_info;
        ++$record_num;

    } #end of while loop

    my $sizeGenes=scalar @mutations;
    printf "The total number of genes extracted from the file %s is %d.\n",$inmaf, $sizeGenes;
    print "\nFor file [$inmaf] we found $record_num mutations\n\n\n";

    close $MAF
        or carp "$0: failed to close infile '$inmaf': $!\n";

    return (\@mutations, \@geneids);
} 
#############################################################################
###################################################################################
#############################################################################
###################################################################################
#############################################################################
###################################################################################
sub Create_pbs_batch
{
    my ($infolder,$outfolder,$scale,$MutRate, $TumorNumber, $TumorType, $gene,$mut_file, $recur_file)=@_;
    my $qsub="";
    my $n= 'pbs_cMACPRF_GeneID_TumorType_';
    print "The format of batch file to be created - bulldogjl (24h/3d) or bulldogk (7d)? \n";
    #my $type= $ARGV[0]; # bulldogjl or bulldogk
    #my $type='bulldogl';
    my $type='bulldogn';

    print "Input the Tumor type: $TumorType\n";
    print "Input the GeneID: $gene\n";

    my $walltime='24h'; # 24h, 3d or 7d

    if ($type eq 'bulldogl' and $walltime eq '24h') { $type='24h_'.$type; }

    elsif ($type eq 'bulldogjl' and $walltime eq '3d') { $type='3d_'.$type; }

    elsif ($type eq 'bulldogk') { $type='7d_'.$type;}
    elsif ($type eq 'bulldogn') { $type=$type;}
    else { print "Error! Input bulldogjl or bulldogk\n";}
    $n=$n.$type.".txt";
    my $infile=$infolder.$n;
    open(IN1,$infile) or croak "Failed to open the file to save pbs batch $infile. Check if you have the required input file pbs_cMACPRF_GeneID_TumorType_bulldogn.txt.\n";
    my @c= <IN1>;
    close(IN1);
    my $c= join("",@c);

    $c=~ s/GeneID/$gene/g;
    $n=~ s/GeneID/$gene/g;
    $c=~ s/TumorcMACPRFMutInput.txt/$mut_file/g;
    $c=~ s/TumorcMACPRFRecurrentInput.txt/$recur_file/g;

    $c=~ s/TumorType/$TumorType/g;
    $n=~ s/TumorType/$TumorType/g;
    
    $c=~ s/TumorSampleNumber/$TumorNumber/g;
    $c=~ s/BackgroundMutationRate/$MutRate/g;
    $c=~ s/Scale/$scale/g;  
    
    my $loc1=index($c,"CSI-MAC -DC ");
    #my $loc2=index($c,"cMACPRF.output.txt",$loc1+2);
    my $loc2=index($c,"CSIMAC.output.txt",$loc1+2);
    my $loc3=index($c,"\n",$loc2+2);
    my $cMACPRF="";
    if ($loc1!=-1 and $loc2!=-1 and $loc3!=-1)
    {
     $cMACPRF="./".substr($c,$loc1,$loc3-$loc1).";\n";
    }
    else {
    printf "Error, empty cMACPRF commands with positions of CSI-MAC -DC %d, output %d, and next line %d!\n", $loc1,$loc2,$loc3;
        
    }
    if ($cMACPRF eq "") { croak "Error, no cMACPRF commands extracted! Check subroutine Create_pbs_batch!\n";}

    printf "Tumor type: %s.\n",$TumorType;
    printf "The file name: %s.\n",$n;
    WriteFile($outfolder,$n,$c);    
    $qsub="qsub $n;\t";
    return ($qsub,$cMACPRF);
}
###############################################################################


###############################################################################
sub TumorNumberCount
{
  my ($mutation_ref)=@_;
  my @mutations = @$mutation_ref;
  my @AllSamples;
  MUTATION:
  foreach my $mutation (@mutations) {
    my $sample=$mutation->{Tumor_Sample_Barcode}; #Lifton_LUSC only
    push @AllSamples,$sample; 
    }
  my $c1=scalar(@AllSamples);
  my @noRedun=RemoveRecurrent(@AllSamples);
  my $sampleCount=scalar (@noRedun);
  #print "Samples after removing redundant:\n",join "\n",@noRedun;
  printf "The number of tumor samples before and after removing redundant: %d/%d\n",$c1,$sampleCount;
  return $sampleCount;
}
###############################################################################


###############################################################################
# For a list, rank it and remove recurrent.
sub RemoveRecurrent
{
   my (@list)=@_;
   @list = sort @list;
   my @new;
   my $t;
   my $size=scalar @list;
   if ($size<2) {@new=@list;}
   for (my $i=1;$i<scalar(@list);$i++)
     {
      $t=$list[$i-1];
      if ($i==1) { push @new, $t;}
      next if $list[$i] eq $t;
      $t=$list[$i];
      push @new, $t; 
     }
   #print "The list before and after removing recurrent:::", scalar @list, "\t", scalar @new, "\n";
   return @new; 
}
###############################################################################

###############################################################################
# For a list, rank it and remove recurrent.
sub KeepRecurrent
{
   my (@list)=@_;
   @list = sort @list;
   my @new;
   my $t;
   my @recurrent;
   my %redun;
   my $count=1;
   for (my $i=0;$i<scalar(@list);$i++)
   {  
      if ($i==0) { push @new, $list[0];}
      else
        {
          #Record the previous one
          $t=$list[$i-1];
          #If the current one is the same as the previous one, Save the recurrent and count
          if ($list[$i] eq $t)
            {
          
              push @recurrent, $t;
              $count=$count+1;
              #printf "Recurrent site: %d\tCount: %d\n",$t,$count;
              #For the last item and is recurrent
              if ($i == scalar(@list)-1)
                {
                  $redun{$t}=$count;
                  if ($new[-1] ne $t) {push @new, $t;}
                }  
              next;
            }
          #if not recurrent, Save the previous recurrent, and empty it for the next   
          else 
            { 
            if ($count>1)
              {
               $redun{$t}=$count;
              }
            @recurrent="";
            $count=1;
            #for the case only two items 
            if ($t ne $list[0]) {
              push @new, $t;          
              }
            #keep the last one if not recurrent
            if ($i==scalar (@list)-1) 
              {
              push @new,$list[$i]; 
              }
            }  
        }
   }
   print "******The list before and after removing recurrent:::", scalar @list, "\t", scalar @new, "\n";
   #print "Before: ",join "\t",@list,"\n";
   #print "After: ",join "\t",@new,"\n";
   print Dumper (\%redun),"\n";
   return (\%redun); 
}

###############################################################################


###############################################################################
#Write contents in the file with the defined file name
sub WriteFile
{
    my ($outfolder,$fileName,$content)=@_;
    $fileName=">".$outfolder.$fileName;
    open (OUT, $fileName);
    my $pres= select(OUT);
    print $content;
    select $pres;
    close(OUT);
}
###############################################################################

###############################################################################
sub openf 
{
    my ($name) = @_;
    open( IN2, $name ) || die "Could not open file $name";
    printf "The input file name for openf: %s.\n", $name;
    my @c = <IN2>;
    close(IN2);
    my $con = join( "", @c );
    #my @l = split( / /, $con );
    #print join ",",@l;
    #return @l;
    return $con;
}
###############################################################################

###############################################################################
sub get_all_isoforms_refflat {
    my ($geneid, $refseq_ref) = @_;
    #Getting all isoforms for the geneid
    my @RefSeq=@$refseq_ref;

    #print "Getting all exons for the gene [$geneid]\n";
    my $gene_rx = qr/$geneid/xm;

    my %coordinates_for;
    #tie %coordinates_for, 'Tie::IxHash';
    my %isoform_length_for;
    my @isoforms;
    my $exon_separator = qq{,};    
    my @strands;
    FEATURE:
    foreach my $datum (@RefSeq) {
        
        my $gene= $datum->{geneName};
        next if ($gene !~ $gene_rx or $gene_rx !~ $gene);

        #record the accession number for each isoform        
        my $accession= $datum->{name};
        push @isoforms, $accession;
        #print "The accession number is $accession.\n";
        
        #Record strand information
        my $strand=$datum->{strand};
        $strand=~ s/ //g;
        if ($strand ne "-" and $strand ne "+") { croak "The strand error $strand should be either + or -\n";}
        push @strands,$strand;
        
        
        my $exon_count=$datum->{exonCount};
        my $TS=$datum->{txStart};
        my $TE=$datum->{txEnd};
        my $cdsS=$datum->{cdsStart};
        my $cdsE=$datum->{cdsEnd};
        my $exonStarts=$datum->{exonStarts};
        my $exonEnds=$datum->{exonEnds};

        my @exonS= split $exon_separator,$exonStarts;
        my @exonE= split $exon_separator,$exonEnds;
        #printf "Isoform %s: CDS start %d; CDS end %d\n",$accession,$cdsS,$cdsE;

        for (my $i=0;$i<$exon_count;$i++)
        {
         my $tmp_s=0;
         my $tmp_e=0;
         #Skip the exon, not start till cdsS
         if ($exonS[$i]<$cdsS and $exonE[$i]<$cdsS) {next; }
         
         #cdsStart and cdsEnd both within one exon
         elsif ($exonS[$i]<=$cdsS and $exonE[$i]>=$cdsS and $exonS[$i]<=$cdsE and $exonE[$i]>=$cdsE)
         {
          $tmp_s=$cdsS;
          $tmp_e=$cdsE;
          $coordinates_for{$accession}{$cdsS} =$cdsE;
          #print "$i+1\t$tmp_s\t$tmp_e\n";
         }  
         
         #cdsStart and exonEnd
         elsif ($exonS[$i]<=$cdsS and $exonE[$i]>=$cdsS and $exonS[$i]<=$cdsE and $exonE[$i]<=$cdsE)
         {
          $tmp_s=$cdsS;
          $tmp_e=$exonE[$i];
          $coordinates_for{$accession}{$cdsS} =$exonE[$i];
          #printf "***Start at CDS! Exon %d: CDS start %d, exon Start %d, exon End %d\n",$i,$tmp_s,$exonS[$i],$tmp_e;
          #print "$i+1\t$tmp_s\t$tmp_e\n";
         }
                
         #exonStart and exonEnd
         elsif ($exonS[$i]>=$cdsS and $exonE[$i]>=$cdsS and $exonS[$i]<=$cdsE and $exonE[$i]<=$cdsE)
         {
          $tmp_s=$exonS[$i];
          $tmp_e=$exonE[$i];
          $coordinates_for{$accession}{$exonS[$i]} =$exonE[$i];
          #print "$i+1\t$tmp_s\t$tmp_e\n";
         } 
         #exonStart and cdsEnd        
          elsif ($exonS[$i]>=$cdsS and $exonE[$i]>=$cdsS and $exonS[$i]<=$cdsE and $exonE[$i]>=$cdsE)
         {
          $tmp_s=$exonS[$i];
          $tmp_e=$cdsE;
          $coordinates_for{$accession}{$exonS[$i]} =$cdsE;
          #printf "***End at CDS! Exon %d: CDS start %d, CDS end %d, exon Start %d, exon End %d\n",$i,$cdsS,$cdsE,$exonS[$i],$exonE[$i];
          #print "$i+1\t$tmp_s\t$tmp_e\n";
          my $exon_length = abs ($tmp_e - $tmp_s);
          $isoform_length_for{$accession} += $exon_length;
          next;
         }     
          else 
          {
          #printf "Isoform %s ends of coordinates at exon %d with total exons %d\n",$accession,$i,$exon_count;
          #printf "---Excluded exons! Exon %d: CDS start %d, CDS end %d, exon Start %d, exon End %d\n",$i,$cdsS,$cdsE,$exonS[$i],$exonE[$i];
          next; }  
              
        #score key stores the isoform number      
         my $exon_length = abs ($tmp_e - $tmp_s);
         $isoform_length_for{$accession} += $exon_length;
                      
        }   
        #printf "Gene length: %d\n\n\n",$isoform_length_for{$accession};      
         
    }
    my $isoform_size=scalar @isoforms;
    
    #for gene not in the refFlat list
    if ($isoform_size==0) { return 0;}
    
    @isoforms = RemoveRecurrent(@isoforms);
    
    #Record the strand information, and it should be the same for the same gene
    @strands= RemoveRecurrent(@strands);
    my $strand_num= scalar(@strands);
    if ($strand_num>1) { 
        print join "\t",@strands,"\n";
        warn "The strand information is not consistant for the same gene $geneid.\n";
        }
    #longest isoform
    my $longest_iso;
    my $longest_ln = 0;
    my $shortest_iso=$isoforms[0];
    my $shortest_ln=$isoform_length_for{$isoforms[0]};
    #print "The shortest isoform is [$shortest_iso] with [$shortest_ln] bp\n";
    foreach my $isoform (sort keys %isoform_length_for) {
        if ($longest_ln < $isoform_length_for{$isoform}) {
            $longest_ln = $isoform_length_for{$isoform};
            $longest_iso = $isoform; 
        }
        if ($shortest_ln > $isoform_length_for{$isoform}) {
            $shortest_ln = $isoform_length_for{$isoform};
            $shortest_iso = $isoform; 
        }         
        
    }
        #Check overlapping for each isoform
    foreach my $isoform (sort keys %coordinates_for) {

        #print "\nChecking for overlapping exons in the isoform [$isoform]\n";
        ##print "******Isoform [$isoform]******\n";
        my %coordinates_isoform = %{$coordinates_for{$isoform}};
        foreach my $start (sort keys %coordinates_isoform)
        {
         ##printf "Start: %d\tEnd: %d\n",$start,$coordinates_isoform{$start};
        }
    
        #print Dumper($coordinates_for{$isoform});
    }
    #print Dumper(\%isoform_length_for);    
    #print "The longest isoform is [$longest_iso] with [$longest_ln] bp\n";
    #print "The shortest isoform is [$shortest_iso] with [$shortest_ln] bp\n";
    return ($longest_ln,$strands[0],\@isoforms,\%coordinates_for,\%isoform_length_for);
}
###############################################################################


###############################################################################
#Extracting all mutations for the given geneid from MAF files, for synonymous and nonsynonymous mutations and keep and remove recurrent positions
sub ExtractMutPos
{
    my ($strand,$ref_coordinates_for,$ref_isoform_length_for,$geneid,$mutation_ref)=@_; 
    my @mutations = @$mutation_ref;

    #we will only plot silent mutations and missense mutations
    #this rx will match Silent and Silent_mutations, same for Missense
    my $interesting_muts_rx = qr/Silent|Missense|Nonsense|Frame_Shift|In_Frame/xmi;
    my $silent_rx = qr/Silent/xmi;
    my $missense_rx = qr/Missense/xmi;
    my $virant_classification_rx = qr/SNP|INS|DEL/xmi;
    my $damaging= qr/Nonsense|Frame_Shift|In_Frame/xmi;
    my %mut_counts;
    my @silent_mut_pos;
    my @missense_mut_pos;
    my @missense_mut_chrom_pos;
    my @damaging_mut_pos;

    #looping through mutations to find different types of mutations for the given geneid
    # to print out missense mutations - the chromosome and the position
    MUTATION:
    foreach my $mutation (@mutations) {
        next MUTATION if $mutation->{Hugo_Symbol} ne $geneid;
    
        if ($mutation->{Variant_Classification} !~ $missense_rx)
        {
            #print "Excluding SNP as Missense Variant: $mutation->{Hugo_Symbol}  $mutation->{Variant_Type} $mutation->{Variant_Classification}\n";
        }
  		printf "variant clasification: %s, variant_type: %s^^^^^^^\n", $virant_classification_rx, $mutation->{Variant_Type};
        #next MUTATION if $mutation->{Variant_Type} !~ $virant_classification_rx;
        #++$mut_counts{SNP};
        #print "Variant: $mutation->{Hugo_Symbol}  $mutation->{Variant_Type} $mutation->{Variant_Classification}\n";
 
        next MUTATION if $mutation->{Variant_Classification} !~ $interesting_muts_rx;
    
        if ($mutation->{Variant_Classification} =~ $silent_rx) {
            push @silent_mut_pos, $mutation->{start};
            my $start= $mutation->{start};
            my $chrom= "Chromosome\t".$mutation->{chr}; 
            my $tmp_strand=$mutation->{strand};
            if ($tmp_strand eq "-") {$tmp_strand=-1;}
            elsif ($tmp_strand eq "+") {$tmp_strand=1;}
            my $ref_allele=$mutation->{ref_allele};
            my $tumor_allele1=$mutation->{newbase};
            my $tumor_allele2=$mutation->{newbase};
            my $tumor_allele=$tumor_allele1;
            if ($tumor_allele1 eq $ref_allele) { $tumor_allele=$tumor_allele2}
            my $tmp_mut_syn=$chrom."\t".$start."\t".$ref_allele."\t".$tumor_allele."\t".$tmp_strand;            
            #print $tmp_mut_syn,"\n";
            
        }
        elsif ($mutation->{Variant_Classification} =~ $missense_rx) {
            push @missense_mut_pos, $mutation->{start};
            my $start= $mutation->{start};
            my $chrom= "Chromosome\t".$mutation->{chr}; 
            my $tmp_strand2=$mutation->{strand};
            if ($tmp_strand2 eq "-") {$tmp_strand2=-1;}
            elsif ($tmp_strand2 eq "+") {$tmp_strand2=1;}
            my $ref_allele=$mutation->{ref_allele};
            my $tumor_allele1=$mutation->{newbase};
            my $tumor_allele2=$mutation->{newbase};
            my $tumor_allele=$tumor_allele1;
            if ($tumor_allele1 eq $ref_allele) { $tumor_allele=$tumor_allele2}
            my $tmp_mut=$chrom."\t".$start."\t".$ref_allele."\t".$tumor_allele."\t".$tmp_strand2;           
            push @missense_mut_chrom_pos, $tmp_mut;
            print $tmp_mut,"\n";
            printf "Chrom: %s\t Position: %s\n",$chrom, $start;
        }

            elsif ($mutation->{Variant_Classification} =~ $damaging) {
            push @damaging_mut_pos, $mutation->{start};
            my $start= $mutation->{start};
            my $chrom= "Chromosome\t".$mutation->{chr}; 
            my $tmp_strand2=$mutation->{strand};
            if ($tmp_strand2 eq "-") {$tmp_strand2=-1;}
            elsif ($tmp_strand2 eq "+") {$tmp_strand2=1;}
            my $ref_allele=$mutation->{ref_allele};
            my $tumor_allele1=$mutation->{newbase};
            my $tumor_allele2=$mutation->{newbase};
            my $tumor_allele=$tumor_allele1;
            if ($tumor_allele1 eq $ref_allele) { $tumor_allele=$tumor_allele2}
            my $tmp_mut=$chrom."\t".$start."\t".$ref_allele."\t".$tumor_allele."\t".$tmp_strand2;           
            #push @damaging_mut_chrom_pos, $tmp_mut;
            #print $tmp_mut,"\n";
            printf "Chrom: %s\t Position: %s\n",$chrom, $start;
        }

        else {
            croak "ERROR\n";
        }

        ++$mut_counts{$mutation->{Variant_Classification}};
    }

    my @noRedun_silent_mut_pos= RemoveRecurrent(@silent_mut_pos);
    my @noRedun_missense_mut_pos= RemoveRecurrent(@missense_mut_pos);
    my @noRedun_missense_mut_chrom_pos=RemoveRecurrent(@missense_mut_chrom_pos);
    my @noRedun_damaging_mut_pos= RemoveRecurrent(@damaging_mut_pos);

    #print "\nTotal mutations for the tumor type for all genes ::: ", scalar @mutations,  "\n";
    print "Total silent mutations for gene $geneid before and after removing recurrents::: ", scalar @silent_mut_pos,"\t", scalar @noRedun_silent_mut_pos,"\n";
    print "Total missense mutations for gene $geneid before and after removing recurrents::: ", scalar @missense_mut_pos,"\t", scalar @noRedun_missense_mut_pos, "\n\n";
    #print Dumper(\%mut_counts);
    print "Recorded missense mutations for gene $geneid before and after removing recurrents::: ", scalar @missense_mut_chrom_pos,"\t", scalar @noRedun_missense_mut_chrom_pos, "\n\n";
    #print join "\n",@noRedun_missense_mut_chrom_pos, "\n";
    print "Total damaging mutations for gene $geneid before and after removing recurrents::: ", scalar @damaging_mut_pos,"\t", scalar @noRedun_damaging_mut_pos, "\n\n";
    
        
    #Convert chromosomal mutation positions to gene CDS positions
    printf "Convert tumor synonymous mutation positions: \n";
    my ($longest_isoform_syn,$max_isoform_syn,$ref_trans_tumor_syn)=TransformMutationPositions($strand, \@silent_mut_pos,$ref_coordinates_for,$ref_isoform_length_for,"S");
    my %trans_tumor_syn=%$ref_trans_tumor_syn;
    #printf "All tumor synonymous mutation positions: \n";
    #print Dumper (\%trans_tumor_syn);
        
    printf "Convert tumor nonsynonymous mutation positions: \n";
    my ($longest_isoform_nonsyn,$max_isoform_nonsyn,$ref_trans_tumor_nonsyn)=TransformMutationPositions($strand,\@missense_mut_pos,$ref_coordinates_for,$ref_isoform_length_for,"R");
    my %trans_tumor_nonsyn=%$ref_trans_tumor_nonsyn;
    printf "All tumor nonsynonymous mutation positions: \n";
    print Dumper (\%$ref_trans_tumor_nonsyn);    


    printf "Convert tumor damaging mutation positions: \n";
    my ($longest_isoform_damaging,$max_isoform_damaging,$ref_trans_tumor_damaging)=TransformMutationPositions($strand,\@damaging_mut_pos,$ref_coordinates_for,$ref_isoform_length_for,"D");
    my %trans_tumor_damaging=%$ref_trans_tumor_damaging;
    printf "All tumor damaging mutation positions: \n";
    print Dumper (\%$ref_trans_tumor_damaging); 



    return ($longest_isoform_syn,$max_isoform_syn,$ref_trans_tumor_syn,$longest_isoform_nonsyn,$max_isoform_nonsyn,$ref_trans_tumor_nonsyn,$longest_isoform_damaging,$max_isoform_damaging,$ref_trans_tumor_damaging);
}
###############################################################################


###############################################################################
# Return the transformed mutation positions for each mutation type (R, S; tumor) from genomic location to gene location for each isoform
# Attention: For negative strand, should start from the end to start or use the total length subtracted
# Include the recurrent sites, remove and keep recurrent later.
sub TransformMutationPositions {
    my ($strand,$ref_mut_pos, $ref_coordinates_for,$ref_isoform_length_for,$MutationType)=@_;
        
    my @mut_pos = @$ref_mut_pos;
    my %coord_for = %$ref_coordinates_for;
    my %isoform_length_for=%$ref_isoform_length_for;
    my %codon_mut_pos;

    @mut_pos = sort {$a<=>$b} @mut_pos;
    #print "\nMutations to be mapped ", scalar @mut_pos, " [@mut_pos]\n";

    # change all isoforms and mutation positions within the corresponding CDS to the new positions as the final gene mutation positions and length
    # keep and return the codon length and position list with the most number of mutations      
    # Go through each isoform
    foreach my $isoform (sort keys %coord_for) {
     
        my %coord_isoform =%{$coord_for{$isoform}};
        #tie %coord_isoform, 'Tie::IxHash';
        my $exonNum=scalar(keys %coord_isoform);
        printf "***Isoform: %s; Gene length: %d; Total exons: %d\n",$isoform,$isoform_length_for{$isoform},$exonNum;
        my $total_mut_map = 0;
        my %mut_seen;
        my %mut_notseen;
        my $exon_num = 0;
        my $new_start = 0;
        my $aa_total_mut_map=0;
        my @length;
        push @length,0;
        my @tmp_mut_pos;
        
        
        # Go through each exons
        # Have to sort keys, to make sure the start positions go from small to big
        foreach my $start (sort keys %coord_isoform)
            {
            $exon_num++;
            my $end=$coord_isoform{$start};       
            my $length = $end-$start; # plus one or not 
            $new_start += $length;
            push @length,$new_start;
            #print "ExonNum: $exon_num\tStart: $start\tEnd: $end\tLength: $length\tNewStart: $new_start\n";
            MUTATION:
            foreach my $mutation (@mut_pos) {
                #recurrent will be removed here.
                #Keep the recurrent sites
                if (exists $mut_seen{$mutation}) {
                    #print "mut $mutation has been seen before\n";
                    ###next MUTATION;
                    }
                if ($mutation <= $end && $mutation >= $start) {
                    ++$aa_total_mut_map;
                    ++$mut_seen{$mutation};
                    my $codon_mut = 0;
                    $codon_mut = $mutation - $start+$length[$exon_num-1];
                    #For negative strand, the position equals total minus assumed positive strand
                    if ($strand eq "-") { 
                        #printf "The strand is negative %s, pos %d.\n",$strand,$codon_mut;
                        $codon_mut=$isoform_length_for{$isoform}-$codon_mut+1;}
                    if ($codon_mut>0 and $codon_mut<$isoform_length_for{$isoform}) {push @tmp_mut_pos,$codon_mut;} # keep only mutations inside the gene isoform
                    ##printf "GenomicMutationPos:\t%d\t GeneMutationPos:\t%d\t ExonStart:\t%d\t ExonEnd:\t%d\t ExonLength:\t%d\t AccumulatedExonLength:\t%d\n",$mutation,$codon_mut,$start,$end,$length,$length[$exon_num-1];
                    } 
                }    
            }
            #printf "%d/", scalar @tmp_mut_pos;
            #my @rm_tmp_mut_pos=RemoveRecurrent(@tmp_mut_pos);
            #printf "%d mutations before/after removing recurrent; \t", scalar @rm_tmp_mut_pos;     
            @{$codon_mut_pos{$isoform}}= @tmp_mut_pos;                  
            foreach my $mutation (@mut_pos) {
                next if exists $mut_seen{$mutation};
                ++$mut_notseen{$mutation};
                #print "Mutation $mutation out of borders\n";
                }           
            my $count_notSeen=scalar(keys %mut_notseen);
            #printf "%d mutations not seen in the isoform\n", $count_notSeen;
        }  

# Get the longest isoform with the most mutations       
       
    #Get the longest isoform
    my $longest_isoform="";
    my $longest_length=0;
    foreach my $isoform (keys %isoform_length_for) {
        my $tmp_length=$isoform_length_for{$isoform};
        if ($tmp_length>$longest_length)
            {
             $longest_length=$tmp_length;
             $longest_isoform=$isoform;
            }
        }
    my $max_isoform=$longest_isoform;
    my @longest_mut_pos=@{$codon_mut_pos{$longest_isoform}};    
    my $longest_mut_size=scalar(@longest_mut_pos);
    my $max_mut=$longest_mut_size; 
            
    foreach my $isoform (sort keys %codon_mut_pos) {
            my @tmp_mut_pos=@{$codon_mut_pos{$isoform}};        
            my $size=scalar(@tmp_mut_pos);
            
            #Get the most mutations
            if ($size>$max_mut) 
                { 
                 $max_mut=$size; 
                 $max_isoform=$isoform; 
                }
            #Get the longest isoform if the number of mutation is the same  
            elsif ($size==$max_mut and $isoform_length_for{$isoform}>$isoform_length_for{$max_isoform}) 
                { 
                 $max_mut=$size; 
                 $max_isoform=$isoform; 
                }
            }
                
                     
    if ($max_isoform eq "") { print "No maximum mutated isoform *** $max_isoform *** found!\n";}       

    printf "The isoform %s with length %d has the maximum mutations %d.\n",$max_isoform,$isoform_length_for{$max_isoform},$max_mut; 
    printf "The longest isoform %s with length %d has the mutations %d.\n",$longest_isoform,$isoform_length_for{$longest_isoform},scalar @{$codon_mut_pos{$longest_isoform}};   
    #printf "Converted Mutation positions: ";
    #print join "\t", @{$codon_mut_pos{$max_isoform}},"\n";  
    return ($longest_isoform,$max_isoform,\%codon_mut_pos); # return codon length as well
}
###############################################################################

#####################################
# Parsing the gene rate output file from MutSigCV
# and returns a ref of hash of hashes
# Attention: 
#####################################

sub MutsigCVparse {
    #Only input MAF filename
    my ($inref) = @_;
    
    # Specify field separator...
    my $RECORD_SEPARATOR = qq{	};
    my $FIELD_COUNT      = 16;

    open my $REF, '<', "$inref"
        or croak "$0: failed to open mutsigcv gene rate file '$inref': $!\n";

    #checking labels of infile
    my $first_line = <$REF>;
    chomp $first_line;
    printf "First line: %s***\n", $first_line;
    my @labels = split $RECORD_SEPARATOR, $first_line, $FIELD_COUNT+1;
    my $size=scalar(@labels);
    my $i=0;
        for my $each (@labels)
    {
    $i++;
    $each =~ s/ //g;
    #printf "**Label %d: %s, %s^^^\n",$i,$each,$ref_labels[$i-1];
    }
        
    #printf "Last item:%s****\n",$labels[$size-1];
    #printf "Fourth item: ***%s****\n",$labels[3];
    my $labels_num = scalar @labels;
    printf "Label number: %d; Field Count: %d\n",$labels_num,$FIELD_COUNT;
    croak "Different number of labels for [$inref]\n"
        if $labels_num != $FIELD_COUNT; 
           
    
    #print "\nFIELDS: [@labels]\n";
    #print "Total fields ::: $labels_num\n";

    my @new_labels = Array::Utils::array_minus(@labels, @MutsigcvLabels);
    #printf "The length of new labels, %d, %s\n",scalar(@new_labels), @new_labels;

    croak "@new_labels labels are not defualt lables" if @new_labels;

    
    my @ref;
    my $record_num = 0;
        while (my $record = <$REF>) {
        next if $record !~ /\S/;
        chomp $record;
        #split columns
        my @values = split $RECORD_SEPARATOR, $record, $FIELD_COUNT+1;
        
        #removing leading and trailing white spaces
        s{^\s+|\s+$}{}g foreach @values;
        
        
        my $values_num = scalar @values;
        my %GeneRates;
        @GeneRates{@MutsigcvLabels} = @values;
               

        #Error if we found a different number of columns
        croak "[$values_num] different from expected [$FIELD_COUNT]" 
            if $values_num != $FIELD_COUNT;

        #Error if there is no Hugo symbol for that gene
        croak "No hugo symbol for record num [$record_num]: $record"
            if !$values[0];
    
        #Populating a hash of array of hashes with records. 
        push @ref, \%GeneRates;
        ++$record_num;

    } #end of while loop

    print "\nFor file [$inref] we found $record_num gene transcripts\n";

    close $REF
        or carp "$0: failed to close infile '$inref': $!\n";

    return \@ref;
    

}



#####################################
# Parsing a whole file of UCSC_refseq_refFlat_hg19
# and returns a ref of hash of hashes
# Attention: chromosome number are not in the standard format (chr1_gl000191_random), double check if needs to use that information.
#####################################

sub parseREFSEQ {
    #Only input MAF filename
    my ($inref) = @_;
    
    # Specify field separator...
    my $RECORD_SEPARATOR = qq{	};
    my $exon_separator = qq{,};
    my $FIELD_COUNT      = 11;

    open my $REF, '<', "$inref"
        or croak "$0: failed to open ref infile '$inref': $!\n";

    #checking labels of infile
    my $first_line = <$REF>;
    chomp $first_line;
    printf "First line: %s***\n", $first_line;
    my @labels = split $RECORD_SEPARATOR, $first_line, $FIELD_COUNT+1;
    my $size=scalar(@labels);
    my $i=0;
        for my $each (@labels)
    {
    $i++;
    $each =~ s/ //g;
    #printf "**Label %d: %s, %s^^^\n",$i,$each,$ref_labels[$i-1];
    }
        
    #printf "Last item:%s****\n",$labels[$size-1];
    #printf "Fourth item: ***%s****\n",$labels[3];
    my $labels_num = scalar @labels;
    printf "Label number: %d; Field Count: %d\n",$labels_num,$FIELD_COUNT;
    croak "Different number of labels for [$inref]\n"
        if $labels_num != $FIELD_COUNT; 
           
    
    #print "\nFIELDS: [@labels]\n";
    #print "Total fields ::: $labels_num\n";

    my @new_labels = Array::Utils::array_minus(@labels, @ref_labels);
    #printf "The length of new labels, %d, %s\n",scalar(@new_labels), @new_labels;

    croak "@new_labels labels are not defualt lables" if @new_labels;

    
    my @ref;
    my $record_num = 0;
        while (my $record = <$REF>) {
        next if $record !~ /\S/;
        chomp $record;
        #split columns
        my @values = split $RECORD_SEPARATOR, $record, $FIELD_COUNT+1;
        
        #removing leading and trailing white spaces
        s{^\s+|\s+$}{}g foreach @values;
        
        
        my $values_num = scalar @values;
        my %gff_field_for;
        @gff_field_for{@ref_labels} = @values;
 

        my @exon_starts= split $exon_separator,$gff_field_for{exonStarts};
        my @exon_ends= split $exon_separator,$gff_field_for{exonEnds};
        #printf "First exon start: %d; Last exon End: %d\n",$exon_starts[0],$exon_ends[-1];
        
        
        croak "Error [$record] with transcription start $gff_field_for{txStart} and end $gff_field_for{txEnd}\n"
            if $gff_field_for{txStart} > $gff_field_for{txEnd};
        croak "Error [$record] with CDS start $gff_field_for{cdsStart} and end $gff_field_for{cdsEnd}\n"
            if $gff_field_for{cdsStart} > $gff_field_for{cdsEnd};
        croak "Error [$record] with exon first start $exon_starts[0] and exon last end $exon_ends[-1]\n"
            if $exon_starts[0] > $exon_ends[-1];  
            
                        
        croak "Error [$record] with transcription start $gff_field_for{txStart} and CDS start $gff_field_for{cdsStart}\n"
            if $gff_field_for{txStart} > $gff_field_for{cdsStart};
        croak "Error [$record] with transcription start $gff_field_for{txStart} and exon first start $exon_starts[0]\n"
            if $gff_field_for{txStart} > $exon_starts[0];   
            
        croak "Error [$record] with transcription end $gff_field_for{txEnd} and CDS end $gff_field_for{cdsEnd}\n"
            if $gff_field_for{txEnd} < $gff_field_for{cdsEnd};        
        croak "Error [$record] with transcription end $gff_field_for{txStart} and exon last end $exon_starts[0]\n"
            if $gff_field_for{txEnd} < $exon_ends[-1]; 
                   
               

        #Error if we found a different number of columns
        croak "[$values_num] different from expected [$FIELD_COUNT]" 
            if $values_num != $FIELD_COUNT;

        #Error if there is no Hugo symbol for that gene
        croak "No hugo symbol for record num [$record_num]: $record"
            if !$values[0];
    
        #Populating a hash of array of hashes with records. 
        push @ref, \%gff_field_for;
        ++$record_num;

    } #end of while loop

    print "\nFor file [$inref] we found $record_num gene transcripts\n";

    close $REF
        or carp "$0: failed to close infile '$inref': $!\n";

    return \@ref;
    

}



#####################################
# Parsing a whole file - Yale 18 columns
# and returns a ref of hash of hashes
#####################################

sub parse_MAF3 {
    #Only input MAF filename
    my ($inmaf) = @_;
    
    # Specify field separator...
    my $RECORD_SEPARATOR = qq{\t};
    #my $FIELD_COUNT      = 20;
    my $FIELD_COUNT      = 18;

    open my $MAF, '<', "$inmaf"
        or croak "$0: failed to open infile '$inmaf': $!\n";

    #checking labels of infile
    my $first_line = <$MAF>;
    $first_line =~ s/\r\n/\n/g; #The problem of Windows (CRLF) \r\n issue
    chomp $first_line;
    printf "First line: %s***\n", $first_line;
    my @labels = split $RECORD_SEPARATOR, $first_line, $FIELD_COUNT+1;
    my $size=scalar(@labels);
    my $i=0;
        for my $each (@labels)
    {
    $i++;
    $each =~ s/ //g;
    #printf "**Label %d: %s, %s^^^\n",$i,$each,$default_labels3[$i-1];
    }
        
    printf "Last item:%s****\n",$labels[$size-1];
    printf "Fourth item: ***%s****\n",$labels[3];
    $labels[3] =~ s/Chromosome/Chrom/g;

    my $labels_num = scalar @labels;
    printf "Label number: %d; Field Count: %d\n",$labels_num,$FIELD_COUNT;
    croak "Different number of labels for [$inmaf]\n"
        if $labels_num != $FIELD_COUNT;
    
    #print "\nFIELDS: [@labels]\n";
    print "Total fields ::: $labels_num\n";

    my @new_labels = Array::Utils::array_minus(@labels, @default_labels3);
    printf "The length of new labels, %d, %s\n",scalar(@new_labels), @new_labels;
    print join "*",@labels,"###\n";
    print join "*",@default_labels3,"***\n";

    croak "***\t@new_labels\t*** labels are not default labels" if @new_labels;
    
    #convert two labels consistent     
    #$default_labels3[14] = 'sample_name1';  # tumor 
    #$default_labels3[15] = 'sample_name2';   # matched normal

    my @mutations;
    my $record_num = 0;
    
    while (my $record = <$MAF>) {
        next if $record !~ /\S/;

        $record=~ s/\r\n/\n/g;
        chomp $record;

        #split columns
        my @values = split $RECORD_SEPARATOR, $record, $FIELD_COUNT+1;
        
        #removing leading and trailing white spaces
        s{^\s+|\s+$}{}g foreach @values;

        #Error if there is no Hugo symbol for that gene
        #croak "No hugo symbol for record num [$record_num]: $record"
        #    if !$values[0];
        #Error if there is no Hugo symbol for that gene

		if (!$values[0]) 
		{
		 print "No hugo symbol for record num [$record_num]: $record\n";
		 next;
		}
       
        #renameing chr X and Y
        #$values[4] = 24 if $values[4] =~ /x/i;
        #$values[4] = 25 if $values[4] =~ /y/i;
        
        #printf "NCBI_Build: *%s***\n", $values[2];
        #renameing chr X and Y
        #printf "Chrom: *%s***\n", $values[3];
        
        #renameing chr X and Y
        $values[3] =~ s/chromosome//i;
        $values[3] =~ s/chr//i;
        $values[3] = uc $values[3] ;
        $values[3] = 24 if $values[3] =~ /x/i;
        $values[3] = 25 if $values[3] =~ /y/i;
        $values[3] = 23 if $values[3] =~ /MT/i;
        $values[3] = 23 if $values[3] =~ /Mito/i;
        
        $values[3] = 'X' if $values[3] =~ /24/;
        $values[3] = 'Y' if $values[3] =~ /25/;
        $values[3] = 'Mito' if $values[3] =~ /23/; 
        
        
        #some mutations are named missense_mutations and other just missense
        
        $values[7] =~ s/_Mutation$//xm;
        $values[8] =~ s/SNV/SNP/xm;
        #printf "Mutations: *%s***\n", $values[7];
        #printf "Value: %s\n";
        
        #Tumor sample name - changing to Patient ID
        #my %barcode_value_for = parse_barcode($values[14]);
        #$values[14] = $barcode_value_for{PatientID};
    ###!!! my %barcode_value_for = MAF::parse_barcode($mutation->{Tumor_Sample_Barcode});
    ###!!! my $sample = $barcode_value_for{PatientID}; # We consider patientID = tumor ID  . '-' . $barcode_value_for{SampleID};

        my $values_num = scalar @values;
        my %record_info;
        @record_info{@default_labels3} = @values;

        #Error if we found a different number of columns
        croak "[$values_num] different from expected [$FIELD_COUNT]" 
            if $values_num != $FIELD_COUNT;

        #Error if there is no Hugo symbol for that gene
       # croak "No hugo symbol for record num [$record_num]: $record"             if !$values[0];

    
        #Populating a hash of array of hashes with records. 
        push @mutations, \%record_info;
        ++$record_num;

    } #end of while loop

    print "\nFor file [$inmaf] we found $record_num mutations\n\n\n";

    close $MAF
        or carp "$0: failed to close infile '$inmaf': $!\n";

    return \@mutations;
} #############################################################################


#####################################
# Parsing a whole file - TCGA 20 columns
# and returns a ref of hash of hashes
#####################################

########CHECK parse_MAF33 for

sub parse_MAF33 {
    #Only input MAF filename
    my ($inmaf) = @_;
    
    # Specify field separator...
    my $RECORD_SEPARATOR = qq{\t};
    my $FIELD_COUNT      = 20;
    #my $FIELD_COUNT      = 18;
    
    open my $MAF, '<', "$inmaf"
        or croak "$0: failed to open infile '$inmaf': $!\n";

    #checking labels of infile
    my $first_line = <$MAF>;
    $first_line =~ s/\r\n/\n/g; #The problem of Windows (CRLF) \r\n issue
    chomp $first_line;
    printf "First line: %s***\n", $first_line;
    my @labels = split $RECORD_SEPARATOR, $first_line, $FIELD_COUNT+1;
    my $size=scalar(@labels);
    my $i=0;
        for my $each (@labels)
    {
    $i++;
    $each =~ s/ //g;
    printf "**Label %d: %s, %s^^^\n",$i,$each,$default_labels33[$i-1];
    }
        
    #printf "Last item:%s****\n",$labels[$size-1];
    #printf "Fourth item: ***%s****\n",$labels[3];
    $labels[3] =~ s/Chromosome/Chrom/g;

    my $labels_num = scalar @labels;
    #printf "Label number: %d; Field Count: %d\n",$labels_num,$FIELD_COUNT;
    croak "Different number of labels for [$inmaf]\n"
        if $labels_num != $FIELD_COUNT;
    
    #print "\nFIELDS: [@labels]\n";
    print "Total fields ::: $labels_num\n";

    my @new_labels = Array::Utils::array_minus(@labels, @default_labels33);
    printf "The length of new labels, %d, %s\n",scalar(@new_labels), @new_labels;

    croak "@new_labels labels are not default labels" if @new_labels;
    
    #convert two labels consistent     
    #$default_labels33[14] = 'sample_name1';  # tumor 
    #$default_labels33[15] = 'sample_name2';   # matched normal
    
    my @mutations;
    my $record_num = 0;
    
    ## GAH Edit: Add list called geneIDs
    my @geneids;

    while (my $record = <$MAF>) {
        next if $record !~ /\S/;
        $record=~ s/\r\n/\n/g;
        chomp $record;
        #split columns
        my @values = split $RECORD_SEPARATOR, $record, $FIELD_COUNT+1;
        
        #removing leading and trailing white spaces
        s{^\s+|\s+$}{}g foreach @values;

        #Error if there is no Hugo symbol for that gene
        #croak "No hugo symbol for record num [$record_num]: $record"
        #    if !$values[0];
        #Error if there is no Hugo symbol for that gene

		if (!$values[0]) 
		{
		 print "No hugo symbol for record num [$record_num]: $record\n";
		 next;
		}
        else {
            push @geneids, $values[0];
        }
		#renameing chr X and Y
        #$values[4] = 24 if $values[4] =~ /x/i;
        #$values[4] = 25 if $values[4] =~ /y/i;
        
        #printf "NCBI_Build: *%s***\n", $values[2];
        #renameing chr X and Y
        #printf "Chrom: *%s***\n", $values[3];
        
        #renameing chr X and Y
        $values[3] =~ s/chromosome//i;
        $values[3] =~ s/chr//i;
        $values[3] = uc $values[3] ;
        $values[3] = 24 if $values[3] =~ /x/i;
        $values[3] = 25 if $values[3] =~ /y/i;
        $values[3] = 23 if $values[3] =~ /MT/i;
        $values[3] = 23 if $values[3] =~ /Mito/i;
        
        $values[3] = 'X' if $values[3] =~ /24/;
        $values[3] = 'Y' if $values[3] =~ /25/;
        $values[3] = 'Mito' if $values[3] =~ /23/; 
        
        
        #some mutations are named missense_mutations and other just missense
        
        $values[7] =~ s/_Mutation$//xm;
        #printf "Mutations: *%s***\n", $values[7];
        #printf "Value: %s\n";
        
        #Tumor sample name - changing to Patient ID
        ###my %barcode_value_for = parse_barcode($values[14]);
        ###$values[14] = $barcode_value_for{PatientID};

        my $values_num = scalar @values;
        my %record_info;
        @record_info{@default_labels33} = @values;

        #Error if we found a different number of columns
        croak "[$values_num] different from expected [$FIELD_COUNT]" 
            if $values_num != $FIELD_COUNT;
                
        #Populating a hash of array of hashes with records. 
        push @mutations, \%record_info;
        ++$record_num;

    } #end of while loop

    print "\nFor file [$inmaf] we found $record_num mutations\n\n\n";

    close $MAF
        or carp "$0: failed to close infile '$inmaf': $!\n";

    return (\@mutations,\@geneids);
} 
#############################################################################






#############################################################################
#
# parsing barcode
#
#BCR aliquot barcode for the tumor sample including the two 
#additional fields indicating plate and well position. 
#i.e. TCGA-SiteID-PatientID-SampleID-PortionID-PlateID-CenterID. 
#The full TCGA Aliquot ID.
#e.g., TCGA-02-0114-01A-01W-0206-08

sub parse_barcode {
    my ($barcode) = @_;
    my %barcode_value_for;

    # Specify field separator...
    my $BARCODE_SEPARATOR = q{-};
    my $FIELD_COUNT       = 7;
    
    my @barcode_fields = qw(TCGA SiteID PatientID SampleID
                            PortionID PlateID CenterID
                           );
    my @values = split $BARCODE_SEPARATOR, $barcode, $FIELD_COUNT+1;
    @barcode_value_for{@barcode_fields} = @values;

    #Some checking
    croak "Not a TCGA barcode [$barcode_value_for{TCGA}]" 
                        if $barcode_value_for{TCGA} ne 'TCGA';

    return %barcode_value_for;
}


#############################################################################
# printing MAF labels in file
#
sub print_MAF_labels_in_file {
    my ($FH) = @_;
    local $" = qq{\t};

    print {$FH} "@default_labels\n";

    return ;
}
#############################################################################

##############################################################################
# printing DATA in MAF file
#
sub print_MAF_file {
    my ($mutation, $FH) = @_;
    local $" = qq{\t};
   
    print {$FH} "@{$mutation}{@default_labels}\n";

    return ;
}
#############################################################################

#############################################################################
sub parse_AA_conservation {
    my ($aa_cons_file) = @_;
    my %cons_for;

    # Specify field separator...
    my $RECORD_SEPARATOR = qq{\t};
    my $FIELD_COUNT      = 8;


    open my $CONS, '<', $aa_cons_file
        or croak "$0: failed to open infile '$aa_cons_file': $!\n";

    #my @SavedItems=("# chr","pos1","Ortholog","Diffspecies");    

    while (my $record = <$CONS>) {
        next if $record !~ /\S/;
        next if $record =~ /^\#/;
        chomp $record;
        
        my @codon = split $RECORD_SEPARATOR, $record, $FIELD_COUNT+1;
        
        #checking if it is codon
        #print "\npedd:: ", abs($codon[1]-$codon[2]-$codon[3]), "\n";
        #croak "Problem with codon positions [@codon]\n"
        #    if abs($codon[1]-$codon[2]) != abs($codon[2]-$codon[3]);
        
        #checking if conservation is a number
        croak "No conservation information for [@codon]\n"
            if $codon[6] !~ /^\d+$/;
       
        foreach my $pos (@codon[1..3]) {
            $cons_for{$pos} = $codon[6]; 
            #$cons_for{$pos} = $codon[0]."_44_".$codon[5]."_".$codon[6];  
            #print $cons_for{$pos}, "\n";          
        }
    }
    close $CONS
        or carp "$0: failed to close infile '$aa_cons_file': $!\n";

    return \%cons_for;

} 
#############################################################################

#############################################################################
sub GetGeneNamesFromRefFlat{
  my ($inf)=@_;
  my $in=openf($inf);
  my @GeneLists=split(/\n/,$in);
  return \@GeneLists;
}
#############################################################################

#Get the complete genome gene length file from refFlat hg19, getting the longest isoform CDS from all isoforms; just need to run once. 
#############################################################################
sub GetGeneLengthFromRefFlat{
  my ($outfolder,$ref_gene_names,$ref_RefFlat)=@_;
  my @GeneList=@$ref_gene_names;
  my %genes_length;
  for my $geneid (@GeneList){
      my ($longest_isoform_length,$strand,$accession_ref,$ref_coordinates_for,$ref_isoform_length_for) = get_all_isoforms_refflat($geneid,$ref_RefFlat);
      $genes_length{$geneid}=$longest_isoform_length;  
  }

  my $size = keys %genes_length;
  printf "The total number genes with length extracted: %d\n",$size;
    open (OUTG, ">".$outfolder."RefFlat_hg19_genelength.txt");
    my $pres= select(OUTG);
    print Dumper(\%genes_length);
    select $pres;
    close(OUTG);
  
  return \%genes_length;
}
#############################################################################

#Get the gene length from the ready processed file
#############################################################################
sub parse_gene_length {
    my ($gl_file) = @_;
    my %genes_length;

   open(IN,$gl_file) or croak "Failed to open the genbank file in parse_genelength_protein_gbk subroutine.\n";
    my @c= <IN>;
    close(IN);
    my $GL= join("",@c);
    #printf "File %s test: \n%s\n",$hg_protein_gbk,substr($genbank,0,100);
    my @GeneL=split(/\n/,$GL);
    my $count=scalar(@GeneL)-1;
    for (my $i=1;$i<scalar(@GeneL);$i++)
    {
        my $l1=index($GeneL[$i],"\'");
        my $l2=index($GeneL[$i],"=>",$l1+1);
        my $l3=index($GeneL[$i],","); # change the file, add the last item with ','
        my $l4=index($GeneL[$i],"}");
        my $tmp_gene="";
        my $tmp_length=0;
        if ($l1!=-1 and $l2!=-1 and $l3!=-1)
        {
        $tmp_gene=substr($GeneL[$i],$l1+1,$l2-$l1-3);
        $tmp_gene=~ s/ //g;
        $tmp_gene= uc $tmp_gene;
        $tmp_length=substr($GeneL[$i],$l2+3,$l3-$l2-3);
        $tmp_length=~ s/\'//g;
        $tmp_length=~ s/ //g;
        $genes_length{$tmp_gene}=$tmp_length;
        }
        else {
          if ($i!=$count){
            warn "parse_gene_length: The locations are -1; double check to confirm!\n";
            printf "The current item %d in the file: %s\n",$i,$GeneL[$i];
          }
          else
          {
            printf "Finish reading all the gene length %d\n",$i-1;
            last;
          }
        }       
    }
    
    my $size = keys %genes_length;
    printf "The total number genes with length extracted: %d\n",$size;
     #print Dumper(\%genes_length);
    return \%genes_length;
}
#############################################################################

#############################################################################
# to do with the same gene with different lengths - assume the same gene are next to each other
sub parse_genelength_protein_gbk {
    my ($hg_protein_gbk) = @_;
    my %genes_length;
    my @genes;
    my @lengths;
    
    
    open(IN,$hg_protein_gbk) or croak "Failed to open the genbank file in parse_genelength_protein_gbk subroutine.\n";
    my @c= <IN>;
    close(IN);
    my $genbank= join("",@c);
    #printf "File %s test: \n%s\n",$hg_protein_gbk,substr($genbank,0,100);
    

    my @gbk=split(/LOCUS       /,$genbank);
    my @syn; # gene with the alternative names
    #printf "The number of genes: %d\n",scalar(@gbk);
    for (my $i=1;$i<scalar(@gbk);$i++)
    
        #for (my $i=0;$i<5;$i++)
    {
     #printf "******\n%s******\n",$gbk[$i];
     my $l1=index($gbk[$i],"                     /gene=\"");
     my $l2=index($gbk[$i],"\"",$l1+28);
     my $l11=index($gbk[$i],"gene_synonym=\"");
     my $l22=index($gbk[$i],";",$l11+1);
     my $l21=index($gbk[$i],"\"",$l11+15);
     my $tmp_syn="";
     
     #for case - more than one syn
     if ($l11!=-1 and $l22!=-1 and $l21!=-1 and $l22<$l21)
     {
          $tmp_syn=substr($gbk[$i],$l11+14,$l22-$l11-14);
          $tmp_syn=~s/ //g;
          $tmp_syn=~s/\n//g;
     }
     # for case - only one syn 
     elsif ($l11!=-1 and $l22==-1 and $l21!=-1)
     {
          $tmp_syn=substr($gbk[$i],$l11+14,$l21-$l11-14);
          $tmp_syn=~s/ //g;
          $tmp_syn=~s/\n//g;
     }
     
     
     if ($tmp_syn ne "")
     {
         push @syn, uc $tmp_syn;
         #printf "Syn: %s\n",$tmp_syn;
     }
 
     
     my $l222=index($gbk[$i],";",$l22+1);
     
     # for case - more than two syns
     while ($l222<$l21 and $l22<$l21 and $l222!=-1 and $l22!=-1 and $l21!=-1 and $l11!=-1)
     {
      $tmp_syn=substr($gbk[$i],$l22+1,$l222-$l22-1);
      $tmp_syn=~s/ //g;
      $tmp_syn=~s/\n//g;
           if ($tmp_syn ne "")
     {
         push @syn, uc $tmp_syn;
         #printf "Syn: %s\n",$tmp_syn;
     }
     $tmp_syn="";
      $l22=$l222;
       $l222=index($gbk[$i],";",$l22+1);    
     }
     
     # for case - more than two syns, save the last item
     if ($l22!=-1 and $l222==-1 and $l21!=-1 and $l11!=-1 and $l22<$l21)
     {
      $tmp_syn=substr($gbk[$i],$l22+1,$l21-$l22-1);
      $tmp_syn=~s/ //g;
      $tmp_syn=~s/\n//g;
     
      if ($tmp_syn ne "")
     {
         push @syn, uc $tmp_syn;
         #printf "Syn: %s\n",$tmp_syn;
     }      
     }       
     my $l3=index($gbk[$i]," aa");
     my $l4=rindex($gbk[$i]," ",$l3-3);
          
     my $tmp_gene="";
     my $length=0;
     
     if ($l1!=-1 and $l2!=-1 and $l3!=-1 and $l4!=-1)
      {
       $tmp_gene=substr($gbk[$i],$l1+28,$l2-$l1-28);
       $tmp_gene=uc $tmp_gene;
       $length=substr($gbk[$i],$l4+1,$l3-$l4-1);
       $length=($length+1)*3;
       $tmp_gene=~s/ //g;
       $length=~s/ //g;

       # Get the longest genes - assuming same name genes are next to each other
       my $tmp_recurS=scalar(@genes);
#
       if ( exists $genes[$tmp_recurS-1] and $tmp_gene =~ m /$genes[$tmp_recurS-1]/)
       {
        #printf "The gene %s recur number %d; last item: %s\n",$tmp_gene,$tmp_recurS,$genes[$tmp_recurS-1];
        #printf Dumper @genes;

         if ($length<$lengths[-1]) {$length=$lengths[-1];}
         }
     
     
        else 
       {
        @genes=();
        @lengths=();
       }
       
       
       #printf "N: %d\tGene: %s\tLength: %d\n",$i+1,$tmp_gene,$length;
       if (exists $genes_length{$tmp_gene})
       {
        my $prev=$genes_length{$tmp_gene}; 
       if ($length>$prev)
       {
       $genes_length{$tmp_gene}=$length;
       }
       }
       else {$genes_length{$tmp_gene}=$length;}
       
       if (scalar(@syn)>=1)
       {
       for (my $j=0;$j<scalar(@syn);$j++)
       {
        my $each= $syn[$j];
        if (exists $genes_length{$each})
        {
        my $prev2=$genes_length{$each}; 
       if ($length>$prev2)
       {
       $genes_length{$each}=$length;
       }
       }
              else {$genes_length{$each}=$length;}
        
       }
       }
        push @genes,$tmp_gene;
        push @lengths, $length;       
           
       }  
       else { last;}    
     
    }
    my $size = keys %genes_length;
    printf "The total number genes with length extracted: %d\n",$size;
     #print Dumper(\%genes_length);
    return \%genes_length;

}
#############################################################################

#############################################################################
sub parse_refgene_dinucleotides {
    my ($dinuc_file) = @_;
    my %dinuc_for;

    # Specify field separator...
    my $RECORD_SEPARATOR = qq{\t};
    my $FIELD_COUNT      = 5;

    open my $DINUC, '<', $dinuc_file
        or croak "$0: failed to open infile '$dinuc_file': $!\n";
    
     my $first_line = <$DINUC>;
     chomp $first_line;
    
     my %gene_seen;
     RECORD:
     while (my $record = <$DINUC>) {
        next if $record !~ /\S/;
        next if $record =~ /^\#/;
        chomp $record;
        
        my @data = split $RECORD_SEPARATOR, $record, $FIELD_COUNT+1;
        
        $data[1] = uc $data[1];        
        my $total_length = $data[2]+$data[3]+$data[4];
        #Duplicated info or different isoform?
        
        if (exists $dinuc_for{$data[1]}) {
            ++$gene_seen{$data[1]};
            #croak "Gene $data[1] duplicated\n" if $gene_seen{$data[1]} > 1;
            #keep only the longest "isoform"
            next RECORD 
                if $total_length <= $dinuc_for{$data[1]}{total};
        }
    
        $dinuc_for{$data[1]}{CpG}       = $data[2]; 
        $dinuc_for{$data[1]}{'TpC/GpA'} = $data[3]; 
        $dinuc_for{$data[1]}{other}     = $data[4];
        $dinuc_for{$data[1]}{total}     = $total_length;

    }
    close $DINUC
        or carp "$0: failed to close infile '$dinuc_file': $!\n";

    return \%dinuc_for;


    return ;
}
#############################################################################
#############################################################################
sub parse_gene_names {
    my ($gl_file) = @_;
    my %genes_length;

   open(IN,$gl_file) or croak "Failed to open the genbank file in parse_genelength_protein_gbk subroutine.\n";
    my @c= <IN>;
    close(IN);
    my $GL= join("",@c);
    #printf "File %s test: \n%s\n",$hg_protein_gbk,substr($genbank,0,100);
    my @GeneL=split(/\n/,$GL);
    for (my $i=1;$i<scalar(@GeneL);$i++)
    {
        my $l1=index($GeneL[$i],"\'");
        my $l2=index($GeneL[$i],"=>",$l1+1);
         my $l22=index($GeneL[$i],"\'",$l2+1);
        #my $l3=index($GeneL[$i],"\'",$l22+2); # change the file, add the last item with ','
        my $l3=rindex($GeneL[$i],",");
        my $tmp_gene="";
        my $tmp_length=0;
        if ($l1!=-1 and $l2!=-1 and $l3!=-1)
        {
        $tmp_gene=substr($GeneL[$i],$l1+1,$l2-$l1-3);
        $tmp_gene=~ s/ //g;
        $tmp_gene= uc $tmp_gene;
        $tmp_length=substr($GeneL[$i],$l2+4,$l3-$l2-5);
        $tmp_length=~ s/\'//g;
        $genes_length{$tmp_gene}=$tmp_length;
        }
        else {
        warn "Parse_gene_names: The locations are -1; double check to confirm!\n";
        printf "Left: %s\n",$GeneL[$i];
        last;
        }       
    }
    
    my $size = keys %genes_length;
    printf "The total number genes with length extracted: %d\n",$size;
     #print Dumper(\%genes_length);
    return \%genes_length;
}
#############################################################################
1;

