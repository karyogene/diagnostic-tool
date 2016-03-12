#!/usr/bin/perl -w
#Script to analyze the pileup files to identify mutations in high-coverage targeted sequencing experiments. If a normal file is present, it only reports somatic mutations. It could also accepts a list of positions to which you want interrogate the pileup file.
#- t parameter (compulsatory): Tumor pileup file.
#-n parameter (optional): Corresponding normal pileup file.
# -p parameter (optional): Positions csv file with header, with the first columns being Chrom,Position. The rest of columns are ignored.
#  -o parameter (compulsory): Output file name.
#use: perl MIDAS.pl -t Tumour.pileup -n Normal.pileup  -o output.csv 
#iv2 23-7-2010., modified Nacho Varela 31-7-2012

use strict;

my $ScriptRun = "\n\nTo run this script use the parameter -t to determine the tumour pileup file name and pathway, the parameter -o for the output file name, 
. Optionally it can accepts the parameters -n to determine the normal pileup file name and pathway and -p for a Positions File to restrict the analysis to it.
\nie: perl MIDAS.pl -t Tumour.pileup -n Normal.pileup  -o output.csv \n\n";

my $TumFile;
my $NormFile;
my %Param;
my $outFile;
my $PositionsFile;




#Import parameters into useful variables. If not present, print a error message.
%Param = @ARGV;
if (!exists($Param{-t}) or !exists($Param{-o})  ){
	print "Error with the parameters:$ScriptRun";
   }
else {
	$TumFile = $Param{-t};
	$outFile = $Param{-o};
}

#check if normal file suministrated and if so, assign it to a variable
if ( exists($Param{-n}) )
{
	$NormFile = $Param{-n};
}

#check if position files has been suministrated and if so, assign it to a variable
if (exists($Param{-p} ) )
{
	$PositionsFile = $Param{-p};
}



#if present a normal file , open it and create a hash containing the data using the Chr_Position info as ID
my %Normlines;
if ( defined $NormFile)
{
	open(NORMFILE, "<$NormFile") or die "Can't open $NormFile\n";
	my @Normlines = <NORMFILE>;
	close NORMFILE;
	for my $j ( 0 .. $#Normlines ){
		my @NormParms = split (/\s+/, $Normlines[$j]);
		my $Name;
		
		#store independently indels and substitutions from the same position
		if ( $NormParms[2] eq '*' )   
		{
			$Name = $Name = $NormParms[0]."_".$NormParms[1]."_".$NormParms[2];
		}
		else
		{
			$Name = $NormParms[0]."_".$NormParms[1];
		}
		
				
		$Normlines{$Name} = $Normlines[$j];
	}
}

#create a variable to store the lines that will be analysed. If positions file present, only those lines will be studied, if not, all the lines in the tumor pileup file will be parsed.
my @lines_to_study;

#if Present Positions file, overwrite the lines_to_study array
if ( defined $PositionsFile)
{
	open(POSITIONS,"<$PositionsFile") or die "Can't open $PositionsFile\n";
	my $header = 0;
	while(my $line = <POSITIONS>)
	{
		#ignore the header
		if ( $header == 0)
		{
			$header++;
			next;
		}
		
		#remove weird characters and add the line to the array
		$line =~ s/\s+//g;
		my @line_items = split(/,/,$line);
		my $name = $line_items[0]."_".$line_items[1];
		push(@lines_to_study,$name);
	}
	close POSITIONS;

}



#Parse similarly the tumor file. 
open (TUMFILE, "<$TumFile") or die "Can't open $TumFile\n";
my @Tumlines = <TUMFILE>;
close TUMFILE;
my %Tumlines;
for my $i ( 0 .. $#Tumlines ){
	my @TumParms = split (/\s+/, $Tumlines[$i]);
	my $Name;
	
	#store independently indels and substitutions from the same position
	if ( $TumParms[2] eq '*' )   
	{
		$Name = $TumParms[0]."_".$TumParms[1]."_".$TumParms[2];
	}
	else
	{
		$Name = $TumParms[0]."_".$TumParms[1];
	}
	
	$Tumlines{$Name} = $Tumlines[$i];
	
	#if not position file present, populate the positions array with the tumor positions
	if ( !defined($PositionsFile) )
	{
		push(@lines_to_study,$Name);
	}
	
}	



#Create the header in the output file
open(OUTFILE, ">>$outFile") or die "Can't write into $outFile\n";
print OUTFILE "Chromosome,Position,WtAllele,MutAllele,Ref_reads,Mut_reads,TumorA,TumorC,TumorG,TumorT,NormalA,NormalC,NormalG,NormalT,MapQualityTum,SeqQualityTum,Result\n";

#Go for each position calculating the alleles and writing the results in the output file
foreach my $Pos(@lines_to_study)
{
	#If the tumor pileup file doesn't contain a line with that position, in means that there is no coverage. Print an empty line and go to the next one.
	if ( !exists($Tumlines{$Pos})  ) {
		my @Pos_items = split(/_/,$Pos);
		print OUTFILE "$Pos_items[0],$Pos_items[1],-,-,-,-,0,0,0,0,0,0,0,0,-,-,No coverage\n";
	}
	else
	{
		my $line = $Tumlines{$Pos};
		my @Tumcols = split (/\s+/, $line);
		my (@Normcols, @NormAlleles,$NormCoverage,@SortedNormAlleles );
		my @NormAlleleCounts = ("-","-","-","-");
		$NormCoverage = 0;
			
		
		#treat in different fashion substitutions and indels
		if ( $Tumcols[2] ne "*" )
		{
				
			#if normal file, get the info from that position
			if ( defined($NormFile) && exists($Normlines{$Pos} ))
			{
				
					@Normcols= split (/\s+/,$Normlines{$Pos});
					@NormAlleles = &calculate_alleles($Normcols[4],$Normcols[5],$Normcols[6],$Normcols[2]);
					$NormCoverage = $NormAlleles[0]{"count"} + $NormAlleles[1]{"count"} + $NormAlleles[2]{"count"} + $NormAlleles[3]{"count"};
					@SortedNormAlleles = sort{ $b->{"count"} <=> $a->{"count"}}@NormAlleles;
					
					#if normal replace the counts in the normal allele counts array
					@NormAlleleCounts = ($NormAlleles[0]{"count"},$NormAlleles[1]{"count"}, $NormAlleles[2]{"count"}, $NormAlleles[3]{"count"});
		
			}
						
			#Calculate the statistics for all the alleles using the subroutine
			my @TumAlleles = &calculate_alleles($Tumcols[4],$Tumcols[5],$Tumcols[6],$Tumcols[2]);
			
				
			#Work out the real enough quality coverage in each sample
			my $TumCoverage = $TumAlleles[0]{"count"} + $TumAlleles[1]{"count"} + $TumAlleles[2]{"count"} + $TumAlleles[3]{"count"};
				
			#Sort the array of hashes by the count of each allele
			my @SortedTumAlleles = sort{ $b->{"count"} <=> $a->{"count"}}@TumAlleles;
			
			my $result = 'ND';
			my $MutAlleleIndex;
			my $MutAllele;
			my $TumMapQuality;
			my $TumSeqQuality;
			my $RefAlleleIndex;
				
			#calculate the index corresponding to the reference allele
			$RefAlleleIndex = 0 if ( $Tumcols[2] eq 'A');
			$RefAlleleIndex = 1 if ( $Tumcols[2] eq 'C');
			$RefAlleleIndex = 2 if ( $Tumcols[2] eq 'G');
			$RefAlleleIndex = 3 if ( $Tumcols[2] eq 'T');
			
				
			
			# If the wt allele is the second in the sorted array, it means that the mut allele is the first element in the sorted array
			if ( $SortedTumAlleles[1]->{"allele"} eq $Tumcols[2] ) {
				$MutAlleleIndex = 0 if ($SortedTumAlleles[0]->{"allele"} eq 'A' );
				$MutAlleleIndex = 1 if ($SortedTumAlleles[0]->{"allele"} eq 'C' );
				$MutAlleleIndex = 2 if ($SortedTumAlleles[0]->{"allele"} eq 'G' );
				$MutAlleleIndex = 3 if ($SortedTumAlleles[0]->{"allele"} eq 'T' );
				$TumMapQuality = $SortedTumAlleles[0]{"MapQualitiesAverage"};
				$TumSeqQuality = $SortedTumAlleles[0]{"SeqQualitiesAverage"};
				$MutAllele = $SortedTumAlleles[0]{"allele"};
			}
			else{
				#If not second allele.
				if ( $SortedTumAlleles[1]{"count"} == 0 ) {
					$MutAllele = '-';
				}
				else {
					$MutAllele = $SortedTumAlleles[1]{"allele"};
				}		
				$MutAlleleIndex = 0 if ($SortedTumAlleles[1]->{"allele"} eq 'A' );
				$MutAlleleIndex = 1 if ($SortedTumAlleles[1]->{"allele"} eq 'C' );
				$MutAlleleIndex = 2 if ($SortedTumAlleles[1]->{"allele"} eq 'G' );
				$MutAlleleIndex = 3 if ($SortedTumAlleles[1]->{"allele"} eq 'T' );
				$TumMapQuality = $SortedTumAlleles[1]{"MapQualitiesAverage"};
				$TumSeqQuality = $SortedTumAlleles[1]{"SeqQualitiesAverage"};
			}
					
				
			#If the coverage is <20 in any of the samples, there is not enough coverage.
			if( ($TumCoverage < 20) ){
				$result = 'No Coverage';
			}
			#If there are not at least 5% of reads supporting a second allele in the tumour the result is Tumor WT.
			elsif ( $TumAlleles[$MutAlleleIndex]{"count"} < 5 || $TumAlleles[$MutAlleleIndex]{"count"} < ($TumCoverage * 0.005) ) {
				$result = 'Tumor WT';
			}
			#If there is evidence of  a third allele the result is 'multiple alleles'
			elsif ( $SortedTumAlleles[2]{"count"} > 5 || $SortedTumAlleles[2]{"count"} >= ($TumCoverage * 0.005) ) {
				$result = 'Multiple alleles';
			}	
			else{
				#If there are calls with enough quality in the tumour and not evidence in the normal it is SOMATIC, if there is evidence in the normal is Germline.
				if ( defined($NormFile) && $NormCoverage > 20 )
				{
					if (  $NormAlleles[$MutAlleleIndex]{"count"} > ($NormCoverage * 0.01)  ) 
					{
						$result = 'germline';
					}
					else
					{
						$result = 'SOMATIC';
						if ( !defined $PositionsFile)
						{
							print OUTFILE "$Tumcols[0],$Tumcols[1],$Tumcols[2],$MutAllele,$TumAlleles[$RefAlleleIndex]{count},$TumAlleles[$MutAlleleIndex]{count},$TumAlleles[0]{count},$TumAlleles[1]{count},$TumAlleles[2]{count},$TumAlleles[3]{count},$NormAlleleCounts[0],$NormAlleleCounts[1],$NormAlleleCounts[2],$NormAlleleCounts[3],$TumMapQuality,$TumSeqQuality,$result\n";
						}
					}
				}
				else {
					$result = 'SUBSTITUTION';
					if ( !defined $PositionsFile)
					{
						print OUTFILE "$Tumcols[0],$Tumcols[1],$Tumcols[2],$MutAllele,$TumAlleles[$RefAlleleIndex]{count},$TumAlleles[$MutAlleleIndex]{count},$TumAlleles[0]{count},$TumAlleles[1]{count},$TumAlleles[2]{count},$TumAlleles[3]{count},$NormAlleleCounts[0],$NormAlleleCounts[1],$NormAlleleCounts[2],$NormAlleleCounts[3],$TumMapQuality,$TumSeqQuality,$result\n";
					}
				}
					
						
			}
			if ( defined $PositionsFile)
			{
				print OUTFILE "$Tumcols[0],$Tumcols[1],$Tumcols[2],$MutAllele,$TumAlleles[$RefAlleleIndex]{count},$TumAlleles[$MutAlleleIndex]{count},$TumAlleles[0]{count},$TumAlleles[1]{count},$TumAlleles[2]{count},$TumAlleles[3]{count},$NormAlleleCounts[0],$NormAlleleCounts[1],$NormAlleleCounts[2],$NormAlleleCounts[3],$TumMapQuality,$TumSeqQuality,$result\n";					
			}	
			
			#This line has been included to correct a bug. If the positions are recovered from a positions file, it is possible that the line in the file corresponds to an indel. In that case, process it as an indel.
			if (defined $PositionsFile and exists($Tumlines{$Pos."_*"}) )
			{
				my $line_2 = $Tumlines{$Pos."_*"};
				my @Tumcols_2 = split(/\s+/,$line_2);
				
				#check the normal hash, if that position is described in the normal with and indel, ignore the line.
				if ( defined($NormFile) && exists($Normlines{$Pos} ) )
				{
					next;
				}		
					
				my $wt_allele;
				my $mut_allele;
				my $coverage = $Tumcols_2[7];
				my $other_reads = $Tumcols_2[12];
				my $mut_reads;
						
						
				#if coverage < 20, ignore 
				
				if ( $coverage < 20 )
				{
					next;
				}	
							
							
				if ($Tumcols_2[8] eq '*')
				{
					$mut_reads = $Tumcols_2[11];
					if ( $Tumcols_2[9] =~ m/\+(\w+)/ )
					{
						$wt_allele = '-';
						$mut_allele = $1;
					}
					elsif ( $Tumcols_2[9] =~ m/\-(\w+)/)
					{
						$wt_allele = $1;
						$mut_allele = '-';
					}	
								
				}
				elsif( $Tumcols_2[9] eq '*')
				{
					$mut_reads = $Tumcols[10];
					if ( $Tumcols_2[8] =~ m/\+(\w+)/ )
					{
						$wt_allele = '-';
						$mut_allele = $1;
					}
					elsif ( $Tumcols_2[8] =~ m/\-(\w+)/)
					{
						$wt_allele = $1;
						$mut_allele = '-';
					}	
								

				}
				else
				{
					next;
				}	
					
				#print only if indels is supported for more than 20 % reads and there are not additional alleles
				if ( ( $mut_reads >= 10 && $mut_reads >= ( 0.05*$coverage) ) && $other_reads == 0 )
				{
					my $position = $Tumcols[1] + 1;
					my $normal_reads = $coverage - $mut_reads;
					print OUTFILE "$Tumcols_2[0],$position,$wt_allele,$mut_allele,$normal_reads,$mut_reads,-,-,-,-,-,-,-,-,-,-,INDEL\n";
						
				}
				
			}
			
			
			
			
			
			
		}
	
				
		# if not it corresponds to a position with an indel
		else
		{
			#check the normal hash, if that position is described in the normal with and indel, ignore the line.
			if ( defined($NormFile) && exists($Normlines{$Pos} ) )
			{
				next;
			}		
				
			my $wt_allele;
			my $mut_allele;
			my $coverage = $Tumcols[7];
			my $other_reads = $Tumcols[12];
			my $mut_reads;
					
					
			#if coverage < 20, ignore 
			
			if ( $coverage < 20 )
			{
				next;
			}	
						
						
			if ($Tumcols[8] eq '*')
			{
				$mut_reads = $Tumcols[11];
				if ( $Tumcols[9] =~ m/\+(\w+)/ )
				{
					$wt_allele = '-';
					$mut_allele = $1;
				}
				elsif ( $Tumcols[9] =~ m/\-(\w+)/)
				{
					$wt_allele = $1;
					$mut_allele = '-';
				}	
							
			}
			elsif( $Tumcols[9] eq '*')
			{
				$mut_reads = $Tumcols[10];
				if ( $Tumcols[8] =~ m/\+(\w+)/ )
				{
					$wt_allele = '-';
					$mut_allele = $1;
				}
				elsif ( $Tumcols[8] =~ m/\-(\w+)/)
				{
					$wt_allele = $1;
					$mut_allele = '-';
				}	
							

			}
			else
			{
				next;
			}	
				
			#print only if indels is supported for more than 20 % reads and there are not additional alleles
			if ( ($mut_reads >= ( 0.2*$coverage) ) && $other_reads == 0 )
			{
				my $position = $Tumcols[1] + 1;
				my $normal_reads = $coverage - $mut_reads;
				print OUTFILE "$Tumcols[0],$position,$wt_allele,$mut_allele,$normal_reads,$mut_reads,-,-,-,-,-,-,-,-,-,-,INDEL\n";
					
			}
				
				


		}
	}	
}
		
	
close OUTFILE;
print "Finished\n";



########################################################################################################################
########################################################################################################################
#Subroutine to calculate allele frequences and qualities. Generate hash for each allele with the count of reads, the average seq quality and the average map quality
########################################################################################################################
sub calculate_alleles {
	my @sequence = split(//, $_[0]);
	my @seqQual = split (//, $_[1]);
	my @mapQual = split (//, $_[2]);
	my $wtBase = $_[3];
	my @alleles = ('A','C','G','T');
	my @counts = ( 0 , 0 , 0 , 0);
	my @SeqQualities = ( 0 , 0 , 0 , 0);
	my @MapQualities = ( 0 , 0 , 0 , 0);
	my @corrected_sequence;
	#These qualities are the base quality + 33 to simplify the comparation. 63 means actually really quality 30.
	my $MinSeqQuality = 58;
	my $MinMapQuality = 48;
	
	
	#This part correct the sequence list to content only one element per base. That is corrected by the indels and special characters that could appear in the sequence string
	my $i = 0;
	while ( $i <= $#sequence) {
		if ( $sequence[$i] =~ m/^[\.\,\*\,a,c,g,t,A,C,G,T]$/ ) {
			$sequence[$i] =~ tr/acgt/ACGT/;
			push(@corrected_sequence,$sequence[$i]);
			$i++;
		}	
		elsif ( ($sequence[$i] eq '+') or ($sequence[$i] eq '-') ) {
			my @slice;
			if ( $i > ($#sequence - 5 )){
				@slice = @sequence[$i .. $#sequence];
			}
			else {
				@slice = @sequence[$i .. ($i+5)];
			}	
			my $substring = join ( '', @slice);
			$substring =~ m/\D(\d+)\D+/g;
			my @digits = split (//, $1);
			my $characters = $#digits + $1 + 2 ;	
			$i = $i + $characters ;
		}
		elsif ( $sequence[$i] eq '^' ) {
			$i += 2;
		}
		
		else {
			$i++;
		}	
	}
	
	#This loop identify the index of the wt allele in the arrays
	my $WtIndex = -1;
	for my $j ( 0 .. $#alleles ){
		if ($alleles[$j] eq $wtBase ){
			$WtIndex = $j;
		}
	}
	
	#~ print "My corrected_sequence array has: $#corrected_sequence + 1 elements\n";
	#~ print "My seqQual array has: $#seqQual + 1 elements\n";
	#~ print "My mapQual array has: $#mapQual + 1 elements\n";
	
	#This part populate the different arrays with the data found in the line.
	for my $k ( 0 .. $#corrected_sequence ) {
		#If . or , use the WT allele
		if ( ($corrected_sequence[$k] eq '.') or ($corrected_sequence[$k] eq ',') ) {
		    $WtIndex < 0 || die "ERROR: No wild-type allele available\n";
				$counts[$WtIndex] += 1;
				$SeqQualities[$WtIndex] += ord($seqQual[$k]) - 33 ;
				$MapQualities[$WtIndex] += ord($mapQual[$k]) - 33 ;
			
			
		}	
		elsif ( ($corrected_sequence[$k] eq 'A') && (ord($mapQual[$k]) > $MinMapQuality ) && ( (ord($seqQual[$k]) > $MinSeqQuality ) ) ){
			
				$counts[0] += 1;
				$MapQualities[0] += ord($mapQual[$k]) - 33;
				$SeqQualities[0] += ord($seqQual[$k]) - 33;	
		}
		elsif ( ($corrected_sequence[$k] eq 'C') && (ord($mapQual[$k]) > $MinMapQuality ) && ( (ord($seqQual[$k]) > $MinSeqQuality ) )) {
			
				$counts[1] += 1;
				$MapQualities[1] += ord($mapQual[$k]) - 33;
				$SeqQualities[1] += ord($seqQual[$k]) - 33;	
		}
		elsif ( ($corrected_sequence[$k] eq 'G') && (ord($mapQual[$k]) > $MinMapQuality ) && ( (ord($seqQual[$k]) > $MinSeqQuality ) )) {
			
				$counts[2] += 1;
				$MapQualities[2] += ord($mapQual[$k]) - 33;
				$SeqQualities[2] += ord($seqQual[$k]) - 33;	
		}
		elsif ( ($corrected_sequence[$k] eq 'T') && (ord($mapQual[$k]) > $MinMapQuality ) && ( (ord($seqQual[$k]) > $MinSeqQuality) )) {
			
				$counts[3] += 1;
				$MapQualities[3] += ord($mapQual[$k]) - 33;
				$SeqQualities[3] += ord($seqQual[$k]) - 33;	
		}
		else {
		}	
		
		
	}
	
	#Finally, we create a hash for each allele to return from the subroutine	
	my %A = ( "allele" => $alleles[0] ,
				    "count" => $counts[0] );
	if ( $A{"count"} == 0 ) {
		$A{"MapQualitiesAverage"} = 0;
		$A{"SeqQualitiesAverage"} = 0;
	}
	else {
		$A{"MapQualitiesAverage"} = $MapQualities[0] / $counts[0];
		$A{"SeqQualitiesAverage"} = $SeqQualities[0] / $counts[0];
	}	
		
	my %C = ( "allele" => $alleles[1] ,
				    "count" => $counts[1]);
	if ( $C{"count"} == 0 ) {
		$C{"MapQualitiesAverage"} = 0;
		$C{"SeqQualitiesAverage"} = 0;
	}
	else {
		$C{"MapQualitiesAverage"} = $MapQualities[1] / $counts[1];
		$C{"SeqQualitiesAverage"} = $SeqQualities[1] / $counts[1];
	}				    
				   
	my %G = ( "allele" => $alleles[2] ,
				    "count" => $counts[2]);
	if ( $G{"count"} == 0 ) {
		$G{"MapQualitiesAverage"} = 0;
		$G{"SeqQualitiesAverage"} = 0;
	}
	else {
		$G{"MapQualitiesAverage"} = $MapQualities[2] / $counts[2];
		$G{"SeqQualitiesAverage"} = $SeqQualities[2] / $counts[2];
	}			    
				    
	my %T = ( "allele" => $alleles[3] ,
				    "count" => $counts[3] );
	if ( $T{"count"} == 0 ) {
		$T{"MapQualitiesAverage"} = 0;
		$T{"SeqQualitiesAverage"} = 0;
	}
	else {
		$T{"MapQualitiesAverage"} = $MapQualities[3] / $counts[3];
		$T{"SeqQualitiesAverage"} = $SeqQualities[3] / $counts[3];
	}				    
			
		
	return ( \%A, \%C, \%G, \%T);
}




	
		
	
	
	


		
		
		
