#!/home/czerik/bin/localperl/bin/perl
#/usr/bin/perl


# This scrict convert the MEME-MAST XML result files into BED file.
	# BED file columns:
		# 1. column: chromosome location of motif
		# 2. column: start of the motif
		# 3. column: end of the motif
		# 4. column: Enhancer ID: enhancer position + Motif name (jaspar ID, alt. name)
		# 5. column: motif p-value
		# 6. column: strand information of motif

# used modules:
use XML::Simple;
use Data::Dumper;
use strict;
use warnings;
use diagnostics;


my $i;
my $Input;
my $e;
my $f;
my $q;

for($i=0;$i<=$#ARGV;$i++){
	if ($ARGV[$i] eq "-in") {$Input = $ARGV[$i+1]}		# MAST xml output file
}

# create object
my $xml = new XML::Simple;
my $data = XMLin("$Input" , ForceArray=>'hit');

# read XML file
my $data_2 = $xml->XMLin("$Input" , ForceArray=>['hit']);
my $motif_lenght= $data_2->{motifs}->{motif}->{width};						# motif lenght will be used in te calculation of te motif ends (genomic position)
my $motif_name= $data_2->{motifs}->{motif}->{name}.'_'.$data_2->{motifs}->{motif}->{best_f};		# Jaspa ID of the motif
my $motif_alt= $data_2->{motifs}->{motif}->{name};						# Alernative motife name 
#Create ID list file for motif mask and BED file creating
open TEMP_OUT, ">${motif_name}_temp.lst" or die"Can't create temporary ID list file\n"; #These temp file will be removed at the end of the analysis. It contains the chr?:????-??? IDs of the enhancer region
	if (ref($data->{sequences})  =~ /ARRAY/) {
		foreach my $ID (@{$data->{sequences}}) {
        	print TEMP_OUT "$_\n" for  keys %{$ID->{sequence}};
}
}
close TEMP_OUT;

#Creating BED files
#open OUT_BED, ">hs_motif_$motif_alt.bed" or die"Can't open hs_motif_$motif_alt.bed\n";		
open LISTAFILE,"${motif_name}_temp.lst" or die "Cant open temp.lst\n";
# if $data->{hit} or {seg} can be ARRAYs ref (if one sequence has more then one motif match)
while (<LISTAFILE>) {
	chomp;
	local $" = "\t";
	my $ID = $_;
	my @ID_pos = split /[:-]+/, $ID;
	if (ref($data_2->{sequences}->{sequence}->{"$ID"}->{seg}) =~ /ARRAY/) {
		foreach  $f (@{$data_2->{sequences}->{sequence}->{"$ID"}->{seg}}) {
                if (ref($f->{hit})  =~ /ARRAY/) {
		foreach  $q (@{$f->{hit}}) {
		print $ID_pos[0]."\t";
                print $ID_pos[1]+$q->{pos}-"1"."\t";
                print $ID_pos[1]+$q->{pos}+$motif_lenght-"1"."\t";
                print $ID."_"."$motif_name"."\t";
                print $q->{pvalue}."\t";
                if (( $q->{strand} eq "forward")) {print "+"} else {print "-"};
                print "\n";
		}
		}
		}
	}
	elsif (ref($data_2->{sequences}->{sequence}->{"$ID"}->{seg}->{hit}) =~ /ARRAY/) {
    		foreach my $e (@{$data_2->{sequences}->{sequence}->{"$ID"}->{seg}->{hit}}) {
		print $ID_pos[0]."\t";
		print $ID_pos[1]+$e->{pos}-"1"."\t";
		print $ID_pos[1]+$e->{pos}+$motif_lenght-"1"."\t";
		print $ID."_"."$motif_name"."\t";
		print $e->{pvalue}."\t";
		if (( $e->{strand} eq "forward")) {print "+"} else {print "-"};
		print "\n";
	}
	} else {
		print $ID_pos[0]."\t";
		print $ID_pos[1]+$e->{pos}-"1"."\t";
		print $ID_pos[1]+$e->{pos}+$motif_lenght-"1"."\t";
		print $ID."_"."$motif_name"."\t";;
		print $e->{pvalue}."\t";
		if (( $e->{strand} eq "forward")) {print  "+"} else {print "-"};
		print "\n";
		}		
}
close LISTAFILE;
#close OUT_BED;
unlink("${motif_name}_temp.lst");
exit;
