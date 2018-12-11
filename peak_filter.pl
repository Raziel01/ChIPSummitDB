#!/usr/bin/perl
use List::Util qw( min max sum);
use Statistics::Basic qw(:all);

for($i=0;$i<=$#ARGV;$i++){
	if ($ARGV[$i] eq "-bed") {$infile1 = $ARGV[$i+1]}
	if ($ARGV[$i] eq "-bedg") {$infile2 = $ARGV[$i+1]}
	if ($ARGV[$i] eq "-of") {$outfile = $ARGV[$i+1]}        #output file}
}

 my @strand_list = qw(pos neg);
foreach $strand (@strand_list){
 
open OUT, ">${outfile}_${strand}.bed" or die"Can't open ${outfile}_${strand}\n";

open (FILE_TEM1, "gunzip -c  $infile2 | awk -v strand=$strand '{OFS=\"\t\"; if     (strand==\"pos\" && \$4>=0) print \$0; else if (strand==\"neg\" &&     \$4\~\/\-\/) print  \$1,\$2,\$3,\$4}' | subtractBed -b stdin -a $infile1  | awk  '{OFS=\"\t\"; print \$1,\$2,\$3,\"0\"}' | sed -e 's/-//g' |") or die "A substr file nem nyithat meg" ;
open (FILE_TEM2, "gunzip -c  $infile2 | awk -v strand=$strand '{OFS=\"\t\"; if (strand==\"pos\" && \$4>=0 && NR>1) print \$0; else if (strand==\"neg\" && \$4\~\/\-\/ && NR>1)  print \$1,\$2,\$3,\$4}' |  sed -e 's/-//g'  |") or die "Az bedgraph file nem nyithat meg";
open (MERGE, ">${outfile}_temp1_merged") or die "Az osszefesult file nem nyithat meg";
my $line1 = <FILE_TEM1>;
my $line2 = <FILE_TEM2>;
while (defined($line1) || defined($line2)) {
	if (defined($line1)) {
	print MERGE $line1;
	$line1 = <FILE_TEM1>;
	}
        if (defined($line2)) {
        print MERGE $line2;
        $line2 = <FILE_TEM2>;
	}
}
close FILE_TEM1;
close FILE_TEM2;
close MERGE;

open (FILE, "awk '{OFS=\"\t\";  print \$1,\$2,\$3,\$1\"_\"\$2\"_\"\$3\"_\"NR,\$4,\$5}' $infile1  |  intersectBed  -a ${outfile}_temp1_merged -b stdin  -wo | sort -k1,1 -k2,2n |  awk '{OFS=\"\t\"; print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' |");


my %hash=();
my %data=();
while (<FILE>) {
	chomp;
	@records2 = split "\t", $_;
	$hash{$records2[7]} = 7;
	($key,$value) = ($records2[7],$records2[3]);
	$data{$key} = [] unless exists $data{$key};
	$interval = ($records2[2] - $records2[1] - 1);
	push @{$data{$key}}, map { $value } 0 .. $interval ;
}
close FILE;

foreach my $key ( keys %data ) {
	my $min = min @{$data{$key}};
	my $max = max @{$data{$key}};
	chomp $max;
	my @count = "";
	push @count, @{$data{$key}};
	chomp @count;
	my @min_index = grep $count[$_]  == 0 , 0.. $#count;
	chomp @min_index;
	my $min_pos = median(@min_index);
	my @max_index = grep $count[$_]  == $max , 0.. $#count;
	chomp @max_index;
	my @min_startend = grep $count[$_]  != $min , 0.. $#count;
	my $min_end = min @{min_startend};	# megkeresi az utolso erteket a peak elott ami nagybb mint nulla (vagy mint a minimum)
	my $min_start = max @{min_startend}; #megkeresi az elso erteket a peak utan ami nagybb mint nulla (vagy mint a minimum)
	$edge1 = "";
	$edge2 = "";
	if (grep $_ == ($min_end - 1), @min_index ) {$edge1 = edge;} else {$edge1 = non}; 
	if (grep $_ == ($min_start + 1) , @min_index) {$edge2 = edge;} else {$edge2 = non};
	my $max_pos = median(@max_index);
	my $median = median(@{$data{$key}});

#Kummulalt osszeg minden sorra percentilis szamitashoz:
$sum=0; 
my @newarr = map { $sum += $_ } @count;

	$sum = eval join '+', @count;
	$firstpercentile_new=($sum/100*25); 
        $secpercentile_new=($sum/100*50);
        $thirdpercentile_new=($sum/100*75);	
my $itr=0;
my $itrb=0;
my $itrc=0;

#Percentile calculating
#
foreach my $number (@newarr){
$itr++ and next if $firstpercentile_new >= $number;
}
foreach my $number (@newarr){
$itrb++ and next if $secpercentile_new >= $number;
}
foreach my $number (@newarr){
$itrc++ and next if $thirdpercentile_new >= $number;
}

# Print out Forward and reverse strand signals
#
$pecarray_1=$newarr[$itr-1];
$pecarray_2=$newarr[$itrb-1];
@key_split = split "_", $key;
$summit_pos=int($key_split[1]+$max_pos);
$filter1_shapefilt = sprintf("%.4f", (($itrc - $itr )/($min_start + 1 - ($min_end - 1))));
$filter1_symmfilt =  sprintf("%.4f",((($itrc -  2 * $itrb + $itr) + 0.001)/(($itrc - $itr)+ 0.001)));
printf OUT  $key_split[0]."\t".$summit_pos."\t".$summit_pos."\t".$key."\t".$filter1_shapefilt."\t".$filter1_symmfilt."\t".$max."\t".$max_pos."\t".$edge1._.$edge2."\n" if exists $hash{$key};

}
%hash=();
`rm -rf ${outfile}_temp1_merged`;
}
close OUT;

#Pairing forward and reverse signals of peaks
#
my %INTER1=();
open INTERFILE1,"${outfile}_pos.bed" or die "Problem with: $Interfilefile1-t\n";
while (<INTERFILE1>){
    s/[\r\t ]+$//g;
    next if (/^$/);
    my $line1 = $_;
    @intersplit1=split('\t',$line1,10);
    $INTER1{$intersplit1[3]} = $line1;
}
close INTERFILE1;

my %INTER2=();
open INTERFILE2,"${outfile}_neg.bed" or die "Problem with: $Interfilefile2-t\n";
while (<INTERFILE2>){
    s/[\r\t ]+$//g;
    next if (/^$/);
    my $line2 = $_;
    @intersplit2=split('\t',$line2,10);
    $INTER2{$intersplit2[3]} = $line2;
}
close INTERFILE2;

open OUT, ">${outfile}_paired.bed" or die"Can't open $outfile\n";

open LISTAFILE,"cat ${outfile}_pos.bed ${outfile}_neg.bed | awk '{print \$4}' | sort | uniq |" or die "Cant open listfile\n";
while (<LISTAFILE>) {
    chomp;
    chomp $INTER2{$_};
    @temp_split2=split('\t',$INTER2{$_});
    chomp $INTER1{$_};
    @temp_split1=split('\t',$INTER1{$_});
    my $focus_ratio = max  (($temp_split2[6]+0.01)/($temp_split1[6]+0.01),($temp_split1[6]+0.01)/($temp_split2[6]+0.01)) ;
    print OUT "$INTER1{$_}\t$INTER2{$_}\t$focus_ratio\n";
}
close LISTAFILE;
close OUT;
`rm -rf   ${outfile}_pos.bed`;
`rm -rf  ${outfile}_neg.bed`; 
exit;


