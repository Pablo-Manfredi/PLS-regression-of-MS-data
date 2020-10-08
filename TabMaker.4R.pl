#!/usr/bin/perl -w
#Name              : TabMaker.4R.pl
#Author            : Dr. Pablo Manfredi
#Goal              : identify SNPs and small INDELS of interest
#Inpout            : all '.AllBase.vcf' files present in the same folder as the script
#Outpout           : a .tsv file containing all candidat alleles emerging during the evolution experiement

use warnings; use Cwd qw(cwd); use strict;

my ($Cuttoff)=1e-6;
my ($Current, @Words, %Phenotype, @headl1, @MSWords, $length, %Conditions, %MSpeptideIntensity);
my ($VcfDIR) = cwd; # get the path of the current directory
opendir (DIR1, $VcfDIR)|| die "can't opendir $VcfDIR: $!";my (@Files_in_Dir) = readdir(DIR1);closedir(DIR1);
foreach my $Files(@Files_in_Dir)
  {
  chomp $Files;
############   PHENOTYPES
  if ($Files =~ /^Phenotype\.tsv$/)
    {
    open(FILEIN, $Files) || die ("Couldn't open the $Files");
    $Current = <FILEIN>;#skip the head line
    Loop1:while (<FILEIN>)
      {
      $Current = $_;
      chomp $Current;
      @Words = split(/\t/, $Current);
      $Phenotype{$Words[0]}=$Words[1];
      }
      close (FILEIN);
    }
############   MS DATA
  elsif ($Files =~ /^MSdata\.tsv$/)
    {
=Example
Pep	        prot      A.1      A.2      A.3      B.1       B.2     B.3
AAFAPQLLDYK	PA2304  6829.45  6076.12  6477.13   9218.43  9464.76  7391.39
AAAAEDPFVISVK	PA3810  8626.04  8873.76  8229.54   9116.45  9343.52  9272.23
GELVQSWPALPAR	PA2867  1439.79  1861.30  1629.01   1295.86  1706.28  1307.05
=cut
    open(FILEIN, $Files) || die ("Couldn't open the $Files");
    $Current = <FILEIN>;#Analyze the head line
    chomp $Current;
    @headl1 = split(/\t/, $Current);
    $length = scalar (@headl1);
    for (my $n =2 ;$n < $length;$n++)
            {
            if ($headl1[$n]  =~/^[^\.]+\.[0-9]+$/)
                {
                $Conditions{$headl1[$n]} = 1;
                }
            }
    while (<FILEIN>)
      {
      $Current = $_;
      chomp $Current;
      @MSWords = split(/\t/, $Current);
      for (my $n =2 ; $n < $length;$n++)
          {
          if ($headl1[$n]  =~/^[^\.]+\.[0-9]+$/)
             {
             $MSpeptideIntensity{$MSWords[1]}{$MSWords[0]}{$headl1[$n]}=$MSWords[$n];
             }
          }
      }
      close (FILEIN);
    }
   }

#################################################################################################
####### Compute Total proteome intensities for normalisation
#################################################################################################
my ($currentProt);
my ($currentAnnot);
my ($currentPep);
my ($condition);
my (%TotalsProteomes);

#############   Compute total MS signal by sample
foreach $currentProt ( keys %MSpeptideIntensity)
  {
  foreach $currentPep ( keys %{$MSpeptideIntensity{$currentProt}})
    {
    foreach $condition ( keys %Conditions)
      {
      if ($MSpeptideIntensity{$currentProt}{$currentPep}{$condition})
        {
        if ($MSpeptideIntensity{$currentProt}{$currentPep}{$condition} =~/^.*NA.*$/){$MSpeptideIntensity{$currentProt}{$currentPep}{$condition} = 0;}
        if ($TotalsProteomes{$condition} && $TotalsProteomes{$condition} != 0){$TotalsProteomes{$condition} = $TotalsProteomes{$condition}+$MSpeptideIntensity{$currentProt}{$currentPep}{$condition};}
        else {$TotalsProteomes{$condition} = $MSpeptideIntensity{$currentProt}{$currentPep}{$condition};}
        }
      }
    }
  }

####### Normalisation to total signal AND give 0.1 TO MISSING VALUES : relative quantification
foreach $condition ( keys %Conditions)
                {
                print $condition.'  '.$TotalsProteomes{$condition}."\n";
                }
foreach $currentProt ( keys %MSpeptideIntensity)
  {
  foreach $currentPep ( keys %{$MSpeptideIntensity{$currentProt}})
    {
    foreach $condition ( keys %Conditions)
      {
      if($MSpeptideIntensity{$currentProt}{$currentPep}{$condition}){$MSpeptideIntensity{$currentProt}{$currentPep}{$condition}=$MSpeptideIntensity{$currentProt}{$currentPep}{$condition}/$TotalsProteomes{$condition};}
      elsif(($MSpeptideIntensity{$currentProt}{$currentPep}{$condition})&&($MSpeptideIntensity{$currentProt}{$currentPep}{$condition} eq 'NA')){$MSpeptideIntensity{$currentProt}{$currentPep}{$condition}=0.1/$TotalsProteomes{$condition};}
      else{$MSpeptideIntensity{$currentProt}{$currentPep}{$condition}=0.1/$TotalsProteomes{$condition};}
      }
    }
  }

#################################################################################################
####### Pooling peptides MS signal per Protein AND get MAXIMA per PROTEIN
#################################################################################################
my (%MSprotIntensity);
foreach $currentProt ( keys %MSpeptideIntensity)
        {
        foreach $currentPep ( keys %{$MSpeptideIntensity{$currentProt}})
                {
                foreach $condition ( keys %Conditions)
                          {
                          if (($MSprotIntensity{$currentProt}{$condition}) && ($MSprotIntensity{$currentProt}{$condition} != 0)){$MSprotIntensity{$currentProt}{$condition} = $MSprotIntensity{$currentProt}{$condition}+$MSpeptideIntensity{$currentProt}{$currentPep}{$condition};}
                          else {$MSprotIntensity{$currentProt}{$condition} = $MSpeptideIntensity{$currentProt}{$currentPep}{$condition};}
                          }
                }
        }

#################################################################################################
####### Preparing Hash for printing
#################################################################################################
my (%Printout);
my (%SUMCHEK);
my ($CurrrrentStrain);
my ($CurrrrentVar);

#Peptides
foreach $currentProt ( keys %MSpeptideIntensity)
        {
        foreach $currentPep ( keys %{$MSpeptideIntensity{$currentProt}})
                {
                foreach $condition ( keys %Conditions)
                          {
                          $Printout{$condition}{$currentPep}=$MSpeptideIntensity{$currentProt}{$currentPep}{$condition};
                          $SUMCHEK{$currentPep}{$condition}=$MSpeptideIntensity{$currentProt}{$currentPep}{$condition};
                          }
                }
        }
#Proteins
foreach $currentProt ( keys %MSprotIntensity)
        {
        foreach $condition (keys %Conditions)
                {
                $Printout{$condition}{$currentProt}=$MSprotIntensity{$currentProt}{$condition};
                $SUMCHEK{$currentProt}{$condition}=$MSprotIntensity{$currentProt}{$condition};
                }
        }
#################################################################################################
####### Comput the maximal value for each variable per file to print (This is used to get rid of variables where all samples are low)
#################################################################################################
my (%MAXCHEKER);
foreach $CurrrrentVar ( keys %SUMCHEK)
        {
        foreach $CurrrrentStrain ( keys %{$SUMCHEK{$CurrrrentVar}})
                        {
                        if ($MAXCHEKER{$CurrrrentVar})
                           {
                           if ($SUMCHEK{$CurrrrentVar}{$CurrrrentStrain}>$MAXCHEKER{$CurrrrentVar})
                              {
                              $MAXCHEKER{$CurrrrentVar} = $SUMCHEK{$CurrrrentVar}{$CurrrrentStrain};
                              }
                           else{}
                           }
                        else
                            {
                            $MAXCHEKER{$CurrrrentVar} = $SUMCHEK{$CurrrrentVar}{$CurrrrentStrain};
                            }
                        }
                }
#################################################################################################
####### Printing out TABLE
#################################################################################################
my ($fileOut) = 'Table4R.tsv';
my (%Counter);

open T4R, '>', $fileOut || die ("Couldn't open the $fileOut file");

#The headline
print T4R "Strain";
print T4R "\t"."Phenotype";
foreach $CurrrrentStrain (sort (keys %Printout))
  {
  foreach $CurrrrentVar (sort (keys %{$Printout{$CurrrrentStrain}}))
    {
    if ($MAXCHEKER{$CurrrrentVar} > $Cuttoff)
      {
      print T4R "\t".$CurrrrentVar;
      }
    }
  }
print T4R "\n";

#The body
my($findBackFirstPart);
my($findBackSecondPart);
my($ValueToPrint);
my($ZerosLength);

foreach $CurrrrentStrain (sort (keys %Printout))
        {
        $Counter{'strains'}++;
        print T4R $CurrrrentStrain;
        if ($Phenotype{$CurrrrentStrain}){print T4R "\t".$Phenotype{$CurrrrentStrain};}else{print T4R "\t".'NA';}
        foreach $CurrrrentVar (sort (keys %{$Printout{$CurrrrentStrain}}))
          {
          if ($MAXCHEKER{$CurrrrentVar} > $Cuttoff)
            {
            $ValueToPrint = sprintf("%.3e", $Printout{$CurrrrentStrain}{$CurrrrentVar});
            if ($ValueToPrint && $ValueToPrint ne '')
              {
              print T4R "\t".$ValueToPrint;
              $Counter{'Values'}++;
              }
            }
          }
         print T4R "\n"; 
         }
close (T4R);
exit;