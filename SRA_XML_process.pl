#!/Users/czerik/local/bin/perl



use XML::Simple;
use Data::Dumper;
#use strict;
#use warnings;
#use diagnostics;


my $i;
my $Input;
my $o;
my $f;
my $q;
my $outfile;
# Input file-ok beolvasasa


if( $ARGV[0] eq '-h' || $ARGV[0] eq '-help')
{
help();
exit;
}


sub help { print '
Useage:

perl /molbio/bin/czerik/SRA_xml_analyser/SRA_XML_process.pl -in input.xml -of output.tbl

	-in	input file, derive from NCBI SRA database. Exported XML file which is a collection of NGS experiments.
	-of	name of output table

Column	HeaderName	Variablename	Description
1.	SRA_ID	$sra_ID	SRA idetifier number
2.	Instrument	$instrument	Sequenator platform
3.	Layout	$layout	Layout: PAIRED/SINGLE
4.	Strategy	$strategy	Method of preparation
5.	Cell_line_and_organism	$cell_line_info	Cell line info: organism abreviation "_" tissue/cell teype "_" cell line
6.	Organism	$organism	Organism
7.	Treatment	$treated	Type of treatment
8.	Treatment_INFO	$treat_tosep	Do the record contain information about treatment?	
9.	Antibody	$antbody	ChIP target
10.	Mutated_factor	$mutation1[0]	Are any factor silenced in cells or not?
11.	Mutation_info	$mutation2	Do the record contain information about silencing?
12.	Mutation_info_source	$mutations_disc	Shows the context, where the mutation are mentioned.
13.	FTPaddress	$ftp_address	Direct download page of raw data.
14.	Warning_messages	@warning_messages	If this column is not empty, the line should be checked manually.

';
}

for($i=0;$i<=$#ARGV;$i++){
            if ($ARGV[$i] eq "-in") {$Input = $ARGV[$i+1]}          # MAST xml output file
        	if ($ARGV[$i] eq "-of") {$outfile = $ARGV[$i+1]}
	}

#Cell line adatbazis megnyitasa. Szukseges a sejttipus azonositasahoz, mert a feltoltok nem egysegesen adtak meg ezt a mezot. Emiatt vagyok kenytelen vagyok keresest hasznalni es nem egyszeruen kinyerni a hash-bol az erteket.        
#open CELL_LINE_FILE, "/Users/czerik/Documents/work/cell_line_project/human_cell_line_database_final3.tbl" or die "Cant open Cell line database\n";
        # create object
open OUT, ">$outfile" or die"Can't open $outfile\n";
my $xml = new XML::Simple;
my $data_2 = $xml->XMLin("$Input" , ForceArray=>['EXPERIMENT_PACKAGE','IDENTIFIERS','PRIMARY_ID']);
print OUT "SRA_ID	Instrument	Layout	Strategy	Cell_line_and_organism	Organism	Treatment	Treatment_INFO	Antibody	Mutated_factor	Mutation_info	Mutation_info_source	FTPaddress	Warning_messages \n";
foreach $o (@{$data_2->{EXPERIMENT_PACKAGE}}) {
    my $sra_ID;
    my $cell_line_info;
    my $treated;
    my $attrubute_num;
    my $antbody;
    my @mutation1;
    my $mutation2;
    my @intersplit2;
    my @antibody;
    my @intersplit3;
    my @mutation1_1;
    my @mutation1_2;
    my @mutation1_3;
    my $mutations_disc;
    my @mutation1_4;
    my $mutation1_5;
    my $j;
    my $value2;
    my $antbody_alt;
    my @antibody_alt;
    my $antibody_alt2;
    my @antibody_alt3;
    my $antbody_alt1;
    my $forj;
    my $aleritem;
    my @warning_messages;
    my $role;
    my %ANTIBODY_HASH;
    my $treat_tosep;
    my %INTER3;
    my @intersplit4;
	my @intersplit5;
    my $antbody_alt_attributenum;
    my %INTER1;
    my @intersplit1;
    my $index_num;
    my $index_cyc;
    my $organism;
    my $organism_ab;
    my @organism;
    my $strategy;
    my $strategy2;
	my @layout;
	my @instrument;
	my $layout;
	my $instrument;
	my $cell_line_string2;
	my $cell_line_string;
	my $ftp_address;
	my $celltype;
	my $cellline;
	my $treated;
	my @real_cell_line;
	my $cellLine_file;
	my $targetprot_file;
	my $attrubute_num7;
    $sra_ID=$o->{EXPERIMENT}->{IDENTIFIERS}->[0]->{PRIMARY_ID}->[0];
    if  ( exists  $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}) { 
        my @cell_line;
#Possible treatment identification
        if (ref($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE})  =~ /ARRAY/) {
            $attrubute_num=@{$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}};

            for $i ( 0 .. $attrubute_num ) {
                push (@cell_line,$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{VALUE});
                    if ($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{TAG}  =~ /treat/i ) {
                        $treated=$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{VALUE};
                        $treated  =~ s/ //g;
                        $treated =~ s/^\s+|\s+$//g;
                        $treated = substr( $treated, 0, 20 );
                    }

                }
		

                if (defined $treated) {
                    $treat_tosep="Yes";
                }
                else {
                    $treated="None";
                    $treat_tosep="No";
                }
            #     print ref($treated);
        }
    elsif (ref($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}) =~ /HASH/) { 
        my $line;
	my  $line2;
        my $elem;
        for $line (keys %{$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}}) {
            push (@cell_line, $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{$line});
    #        push (@cell_line,$elem);
        }   
	if ( exists $o->{EXPERIMENT}->{EXPERIMENT_ATTRIBUTES} && ref($o->{EXPERIMENT}->{EXPERIMENT_ATTRIBUTES}->{EXPERIMENT_ATTRIBUTE})  =~ /HASH/) {
	  for $line2 (keys %{$o->{EXPERIMENT}->{EXPERIMENT_ATTRIBUTES}->{EXPERIMENT_ATTRIBUTE}}) {
            push (@cell_line, $o->{EXPERIMENT}->{EXPERIMENT_ATTRIBUTES}->{EXPERIMENT_ATTRIBUTE}->{$line2});
        }
	}        
            if ($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG}  =~ /treat/i ) {
                $treated=$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{VALUE};
                $treated =~ s/ //g;
                $treated = substr( $treated, 0, 20 );
            }   
            if (defined $treated) {
                print $treated;
            }
            else {
                $treated="None";
            }

    }
	
	if  (exists $o->{EXPERIMENT}->{EXPERIMENT_ATTRIBUTES} && ref($o->{EXPERIMENT}->{EXPERIMENT_ATTRIBUTES}->{EXPERIMENT_ATTRIBUTE})  =~ /ARRAY/) {
		$attrubute_num7=@{$o->{EXPERIMENT}->{EXPERIMENT_ATTRIBUTES}->{EXPERIMENT_ATTRIBUTE}};
		for $i ( 0 .. $attrubute_num7 ) {
		push (@cell_line,$o->{EXPERIMENT}->{EXPERIMENT_ATTRIBUTES}->{EXPERIMENT_ATTRIBUTE}->[$i]->{TAG});
		push (@cell_line,$o->{EXPERIMENT}->{EXPERIMENT_ATTRIBUTES}->{EXPERIMENT_ATTRIBUTE}->[$i]->{VALUE});
		}
	}
    push(@cell_line, ($o->{SAMPLE}->{TITLE}, $o->{SAMPLE}->{alias}, $o->{EXPERIMENT}->{alias}));                              
    $cell_line_string = join(" ", @cell_line);
    $cell_line_string =~ s/_/ /g;
#    $cell_line_string =~ s/-/ /g;
    $cell_line_string =~ s/:/ /g;
    $cell_line_string =~ s/,/ /g;
#    $cell_line_string =~ s/[#\-%&\$*+()]/ /g;
$cell_line_string =~ s/cell //gi;
$cell_line_string =~ s/cells //gi;
$cell_line_string =~ s/cells,//gi;
$cell_line_string =~ s/cell,//gi;
$cell_line_string =~ s/male //gi;
$cell_line_string =~ s/female //gi;
$cell_line_string =~ s/control //gi;
$cell_line_string =~ s/cell//gi;
$cell_line_string =~ s/line//gi;
$cell_line_string =~ s/match//gi;
$cell_line_string =~ s/study//gi;
$cell_line_string =~ s/sample//gi;
$cell_line_string =~ s/from//gi;
$cell_line_string =~ s/immortaliz//gi;
$cell_line_string =~ s/therapy//gi;
$cell_line_string =~ s/elapse//gi;
$cell_line_string =~ s/notapplicable//gi;
$cell_line_string =~ s/France//gi;
$cell_line_string =~ s/treatment//gi;
$cell_line_string2=$cell_line_string;
#$cell_line_string2 =~ s/-//g;
$cell_line_string =~ s/-/ /g;
$cell_line_string =~ s/[#\-%&\$*+()]/ /g;
    #new_if wron delete following line
    $cell_line_string = "$cell_line_string $o->{SAMPLE}->{TITLE} $o->{SAMPLE}->{alias} $o->{EXPERIMENT}->{alias}";
    
#Organism identification    
    $organism=$o->{SAMPLE}->{SAMPLE_NAME}->{SCIENTIFIC_NAME};
	$organism =~ s/ /_/g;
    if ($organism =~ /homo/i && $organism =~ /sapiens/i) {
        $organism_ab="hs";
    }
    elsif ($organism =~ /mus/i && $organism =~ /musculus/i) {
        $organism_ab="mm";
    }
    if ($organism eq "") {
        if ($cell_line_string =~ /human/i || $cell_line_string =~ /homo/i && $cell_line_string =~ /sapiens/i) {
            $organism_ab="hs";
            $organism="Homo sapiens";
        }
        elsif ($cell_line_string =~ /mouse/i || $cell_line_string =~ /mice/i || $cell_line_string =~ /mus/i && $cell_line_string =~ /musculus/i) {
            $organism_ab="mm";
            $organism="Mus Musculus";
        }
        else {
            $organism="Unknown";
            $organism_ab="Unknown";
        }
    }
   #$organism =~ s/ /_/g;
	 elsif ($organism ne "" && $organism_ab ne "hs" && $organism_ab ne "mm"  ) {
            @organism=split("_", $organism,3);
		#@organism=split(/^a-zA-Z\d\s:/,$organism,2);
            $organism_ab=substr( $organism[0], 0, 1 ).substr( $organism[1], 0, 1 );
            push (@warning_messages,'|The results of this organism type is not processable yet, the experiment result are unrliable|');
    }

	if ($organism_ab eq "hs") {
	$cellLine_file="/molbio/bin/czerik/cell_line_db/hs/human_cell_lineDB.tbl";
	$targetprot_file="/molbio/bin/czerik/DNA_assciated_protein_db/hs/DNA_associated_database_final3.tbl";
	}

	elsif ($organism_ab eq "mm") {
        $cellLine_file="/molbio/bin/czerik/cell_line_db/mm/mouse_cell_lineDB.tbl";
        $targetprot_file="/molbio/bin/czerik/DNA_assciated_protein_db/mm/DNA_associated_database_final3.tbl";
        }

	else {
	$cellLine_file="/molbio/bin/czerik/cell_line_db/hs/human_cell_lineDB.tbl";
        $targetprot_file="/molbio/bin/czerik/DNA_assciated_protein_db/hs/DNA_associated_database_final3.tbl";
	}

  # $organism =~ s/ /_/g; 
    #my @real_cell_line;
    if ($cell_line_string =~ /IgG/) {
        push (@warning_messages,'|The Experiment can be IgG control|');
    }
#Strategy identification
    #$strategy=$o->{EXPERIMENT}->{DESIGN}->{LIBRARY_DESCRIPTION}->{LIBRARY_STRATEGY};
#    $strategy=$o->{EXPERIMENT}->{DESIGN}->{LIBRARY_DESCRIPTOR}->{LIBRARY_STRATEGY};
	
	push (@layout, keys %{$o->{EXPERIMENT}->{DESIGN}->{LIBRARY_DESCRIPTOR}->{LIBRARY_LAYOUT}});
	$layout=$layout[0];
	if ($layout eq "") {
		$layout="Unknown";
	}

	push(@instrument,keys %{$o->{EXPERIMENT}->{PLATFORM}});
	$instrument=$instrument[0];
	if ($instrument eq "") {
		$instrument="Unknown";
	}

	

	$strategy=$o->{EXPERIMENT}->{DESIGN}->{LIBRARY_DESCRIPTOR}->{LIBRARY_STRATEGY};

    if ($strategy =~ /chip/i && $strategy =~ /seq/i) {
        $strategy2="ChIP-seq";
    }
    elsif ($strategy =~ /chip/i && $strategy =~ /exo/i) {
        $strategy2="ChIP-exo";
    }
    elsif ($strategy =~ /rna/i && $strategy =~ /seq/i) {
        $strategy2="RNA-seq";
    }
    elsif ($strategy =~ /gro/i && $strategy =~ /seq/i) {
        $strategy2="GRO-seq";
    }
    elsif ($strategy =~ /dnase/i && $strategy =~ /seq/i) {
        $strategy2="DNase-seq";
    }
    elsif ($strategy =~ /ChIA/i && $strategy =~ /pet/i) {
        $strategy2="ChIA-PET";
    }   
    elsif ($strategy =~ /FAIRE/i && $strategy =~ /seq/i) {
        $strategy2="FAIRE-seq";
    }   
    elsif ($strategy =~ /HiC/i && $strategy =~ /seq/i) {
        $strategy2="HiC-seq";
    }
	elsif (defined $strategy && !defined $strategy2)  {
		$strategy2=$strategy;
	}
    if ($strategy eq "" && $cell_line_string =~ /chip-seq/i ||  $cell_line_string =~ /chip/i && $cell_line_string =~ /seq/i ) {
        $strategy="ChIP-seq";
        $strategy2="ChIP-seq";
    }
    elsif  ($strategy eq "" && $cell_line_string =~ /chip-exo/i || $cell_line_string =~ /chip/i && $cell_line_string =~ /exo/i) {
        $strategy="ChIP-exo";
        $strategy2="ChIP-exo";
    }
    elsif  ($strategy eq "" && $cell_line_string =~ /rna/i && $cell_line_string =~ /seq/i) {
        $strategy="RNA-seq";
        $strategy2="RNA-seq";
    }
    elsif  ($strategy eq "" && $cell_line_string =~ /dnase/i && $cell_line_string =~ /seq/i) {
        $strategy="DNase-seq";
        $strategy2="DNase-seq";
    }
    elsif  ($strategy eq "" && $cell_line_string =~ /ChIA/i && $cell_line_string =~ /pet/i) {
        $strategy="ChIA-PET";
        $strategy2="ChIA-PET";
    }
    elsif  ($strategy eq "" && $cell_line_string =~ /FAIRE/i && $cell_line_string =~ /seq/i) {
        $strategy="FAIRE-seq";
        $strategy2="FAIRE-seq";
    }
    elsif  ($strategy eq "" && $cell_line_string =~ /gro/i && $cell_line_string =~ /seq/i) {
        $strategy="RNA-seq";
        $strategy2="RNA-seq";
    }
    elsif  ($strategy eq "" && $cell_line_string =~ /HiC/i && $cell_line_string =~ /seq/i) {
        $strategy="HiC-seq";
        $strategy2="HiC-seq";
    }
	elsif  ($strategy eq "" && $cell_line_string =~ /Hi-C/i ) {
		$strategy="HiC-seq";
		$strategy2="HiC-seq";
	}
    if  (!defined $strategy ) {
        $strategy="Unknown";
        $strategy2="Unknown";
        }



    if ($strategy2 ne "ChIP-seq" && $strategy2 ne "ChIP-exo" && $strategy2 ne "ChIA-PET"){
        $antbody=$strategy2;
    }
    
#Identifying antibody and possible mutation
    open FACTOR_FILE, $targetprot_file or die "Cant open Factor list database\n";
    if ( $o->{SAMPLE}->{TITLE}  =~ /mut/ || $o->{SAMPLE}->{alias} =~ /mut/ ||  $o->{SAMPLE}->{TITLE}  =~ /_kd/ || $o->{SAMPLE}->{alias} =~ /_kd/ || $cell_line_string =~ /si[A-Z]/ || $cell_line_string =~ /sh[A-Z]/ ||  $cell_line_string =~ /KO/ ||  $cell_line_string =~ /knock/i &&  $cell_line_string  =~ /OUT/i ) {
        @intersplit2=split(' ',$cell_line_string,200);
        @intersplit2 = grep defined, @intersplit2;
        @mutation1_1 =grep(/si.[A-Z].*/o,@intersplit2);
        @mutation1_4 =grep(/si[A-Z]./o,@intersplit2);
        @mutation1_2 =grep(/sh.[A-Z].*/o,@intersplit2);
        @mutation1_3 =grep(/[A-Za-z].*._mut_.*/o,@intersplit2);
        #@mutation1_5 =grep(/[A-Z].*._kd.*/o,@intersplit2);
        if ($cell_line_string  =~ /_kd/ || $cell_line_string =~ /knock/i &&  $cell_line_string  =~ /OUT/i ||  $cell_line_string  =~ /down/i ||  $cell_line_string  =~ /KO/)  {
            $index_num=@intersplit2;
            for $index_cyc (1 .. $index_num) {
                #print  $index_cyc;
                if ($intersplit2[$index_cyc] eq "kd" || $intersplit2[$index_cyc] eq "ko") {
                    #     print "$intersplit2[$index_cyc] $index_cyc";
                    if (defined $intersplit2[$index_cyc-1]){
			$mutation1_5=$intersplit2[$index_cyc-1];
                	}
		}
		if ($mutation1_5 eq "") {
		$mutation1_5="NA";
		}
            }
        }
        push( @mutation1, (@mutation1_3, @mutation1_2, @mutation1_1, @mutation1_4,$mutation1_5) );
        $mutation2 = "YES";
        if (! @mutation1){
            push(@mutation1,"NA");
        }
    }
    if (! @mutation1) {
        push (@mutation1,"NA");
    }
    if (!defined $mutation2) {
        $mutation2="NA";
    }
	@mutation1= grep defined, @mutation1;
    $mutations_disc=$mutation1[0];
    $mutation1[0]=~ s/si//g;
    $mutation1[0]=~ s/sh//g;
    $mutation1[0]=~ s/_mut_(\w+)//g;
    $mutation1[0]=~ s/_kd(\w+)//g;
    $mutation1[0] = uc $mutation1[0];

    if ( $o->{SAMPLE}->{TITLE}  =~ /input/i ||  $o->{EXPERIMENT}->{alias}  =~ /input/i ) {
        $antbody="INPUT"; 
    }
#    if ( $o->{SAMPLE}->{TITLE}  =~ /IgG/ ||  $o->{EXPERIMENT}->{alias}  =~ /IgG/ ) {
#        $antbody="IgG";
#    }
    if ( $antbody  eq "") {
    while (<FACTOR_FILE>) {
        if ( $o->{SAMPLE}->{TITLE}  =~ /IgG/i ||  $o->{EXPERIMENT}->{alias}  =~ /IgG/i ) {
		 push (@antibody,"IgG");
	}
	chomp;
        my $line2 = $_;
        @intersplit3=split('\t',$line2,2);
        $INTER3{$intersplit3[0]} = $intersplit3[1];
	if  ($antbody eq "" ) {
		if (ref($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE})  =~ /ARRAY/) {
		for $j ( 0 .. $attrubute_num ) {
			if  ($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$j]->{TAG}  eq "antibody" || $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{TAG}  eq "chip antibody" ||  $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{TAG}  eq "hgn" && $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$j]->{TAG}  !~ /vendor/i) {
				if ($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$j]->{VALUE} =~ /(?:\W|^)(\Q$intersplit3[0]\E)(?:\W|$)/i) {
					push (@antibody,$intersplit3[0]);
				}
			}
		}
		}
	}
	if ( ! @antibody && $cell_line_string2 =~ /anti-/i) {
        @intersplit5=split(' ',$cell_line_string2,200);
        @intersplit5 = grep defined, @intersplit5;	
	@antibody =grep(/anti-.[A-Za-z].*/oi,@intersplit5);
	$antibody[0]=~ s/anti-//gi;
	$antbody=$antibody[0];
	#print join (" ", @intersplit5);
	}
        elsif ( ! @antibody && $cell_line_string =~ /(?:\W|^)(\Q$intersplit3[0]\E)(?:\W|$)/i) {
            push (@antibody,$intersplit3[0]);
        }
        elsif (! @antibody && $cell_line_string =~ /$intersplit3[0]/i) {
                @intersplit4=split(' ',$cell_line_string,200);
                #push(@antibody_alt,$intersplit3[0].(grep(/$intersplit3[0]/i,@intersplit4))) ;            
                @antibody_alt =grep(/$intersplit3[0]/i,@intersplit4); 
                for my $aleritem (@antibody_alt) {
                    $ANTIBODY_HASH{$intersplit3[0]} = $aleritem;
                    #$aleritem2="$aleritem-$intersplit3[0]";
                    push (@antibody_alt3, keys %ANTIBODY_HASH);
                    
                }
            }
            # print join(" ",@antibody_alt3 );
                
    }
    #my $antbody;
    for my $value1 (@antibody) {
        if (length $value1 > length $antbody && $value1 !~ $mutation1[0]) {
            $antbody = $value1;
        }
    }
    $antbody_alt_attributenum=@antibody_alt3;
    for $forj (0 .. $antbody_alt_attributenum) {
        for $role ($antibody_alt3[$forj]) {
        $value2=$role;
        if (length $value2 > length $antbody_alt1 && $value2 !~ $mutation1[0]) {
            $antbody_alt1 = $value2;
        }
    }
    }
    $antbody_alt="$antbody_alt1-$ANTIBODY_HASH{$antbody_alt1}";
    #$antbody_alt="$antbody_alt1-$antibody_alt3->[$forj]->{$antbody_alt1}";
    if  ($antbody eq "" && defined $antbody_alt ) { 
        $antbody=$antbody_alt;
        push (@warning_messages,'|Unusual antibody target|');
    }

    if  ($antbody eq "" ) {
        if (ref($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE})  =~ /ARRAY/) {
            for $j ( 0 .. $attrubute_num ) {
            if  ($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$j]->{TAG}  eq "antibody" || $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{TAG}  eq "chip antibody"  && $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$j]->{TAG}  !~ /vendor/i) {
            $antbody=$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$j]->{VALUE};
            $antbody=~ s/abcam//gi;
            $antbody=~ s/ab[1-9]//g;
            $antbody=~ s#\(.*\)##g;
            $antbody=~ s/Santa|Cruz|Biotechnology//gi;
            $antbody=~ s/^\s+|\s+$//g;
            $antbody=substr( $antbody, 0, 10 );
            }
            }
        }
    }
        elsif (ref($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}) =~ /HASH/) {
            if ($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG}  !~ /vendor/i && $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG}  eq "antibody" || $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG}  =~ /chip/i && $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG}  =~ /antibody/i) {
                $antbody=$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{VALUE};
                $antbody=~ s/abcam//gi;
                $antbody=~ s/ab[1-9]//g;
                $antbody=~ s#\(.*\)##g;
                $antbody=~ s/Santa|Cruz|Biotechnology//gi;
                $antbody=~ s/^\s+|\s+$//g;
                $antbody=substr( $antbody, 0, 10 );
                
            }

    }
    }
    if  ($antbody eq "" && $cell_line_string =~ /input/i ) {
        $antbody="INPUT";
    }
    if  ($antbody eq "" ) {
        $antbody="missing";   
    }
#Identifying origins
    open CELL_LINE_FILE, $cellLine_file or die "Cant open Cell line database\n";
    while (<CELL_LINE_FILE>) {
        chomp;
        my $line1 = $_;
	$cell_line_string2 =~ s/-//g;
	$cell_line_string2 =~ s/[#\-%&\$*+()]/ /g;
        @intersplit1=split('\t',$line1,3);
        $INTER1{$intersplit1[0]} = $intersplit1[1];
        if ($cell_line_string2 =~ /(?:\W|^)(\Q$intersplit1[0]\E)(?:\W|$)/i) {
            push (@real_cell_line,$intersplit1[0]);
        }
    }
    my $cell_line_final;
    for my $value (@real_cell_line) {
        if (length $value > length $cell_line_final) {
            $cell_line_final = $value;
        }
    }
    $cell_line_info=$INTER1{$cell_line_final};

#Ha nem ismert a sejtvonal
#    my $celltype;
#    my $cellline;
#    my $treated;
    unless (defined $cell_line_final) {
        if (ref($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE})  =~ /ARRAY/) {
            $attrubute_num=@{$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}};
            for $i ( 0 .. $attrubute_num ) {
                if ($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{TAG}  =~ /tissue/i || $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{TAG}  =~ /name/i && $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{TAG}  =~ /source/i ||$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{TAG}  =~ /cell/i && $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{TAG}  =~ /type/i ) {
                    $celltype=$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{VALUE};
                    $celltype =~ s/ //g;
                    $celltype = substr( $celltype, 0, 20 );
                }
                elsif ($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{TAG}  =~ /cell/i && $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{TAG}  =~ /line/i  ) {
                    $cellline=$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->[$i]->{VALUE};
                    $cellline =~ s/ //g;
                    $cellline = substr( $cellline, 0, 10 );
                    # print $cellline;
                }
            }
        }
        elsif (ref($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}) =~ /HASH/) {
            if ( $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG}  =~ /tissue/i || $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG}  =~ /name/i && $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG} =~ /source/i || $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG}  =~ /cell/i && $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG} =~ /type/i ) {
                $celltype=$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{VALUE};
                $celltype =~ s/ //g;
                $celltype = substr( $celltype, 0, 20 );
                
            }
            elsif ($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG}  =~ /cell/i && $o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG} =~ /line/i ) {
                $cellline=$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{VALUE};
                $cellline =~ s/ //g;
                $cellline = substr( $celltype, 0, 20 );
            }

            elsif ($o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{TAG}  =~ /treat/i ) {
                $treated=$o->{SAMPLE}->{SAMPLE_ATTRIBUTES}->{SAMPLE_ATTRIBUTE}->{VALUE};
                $treated =~ s/ //g;
                $treated = substr( $treated, 0, 20 );
            }

        }
	$celltype =~ s/_/ /g;   
	$cellline =~ s/_/ /g;
        if ( defined $celltype && defined  $cellline) {
            $cell_line_info="$organism_ab\_$celltype\_$cellline";
        }
        if ($celltype eq "" && $cellline eq "" ) {
            $cell_line_info="$organism_ab\_Unknown_Unknown";
        }
        if  (defined $celltype && !defined  $cellline ) {
            $cell_line_info="$organism_ab\_$celltype\_$celltype";
        }
        if  (!defined $celltype && defined  $cellline ) {
            $cell_line_info="$organism_ab\_$cellline\_$cellline";
        }
        if  ($treated eq "" ) {
            $treated="none";
        }
        } 
    $cell_line_info =~ s/[#\-%&\$*+()]/ /g;
	$cell_line_info =~ s/ //g;
    }

    else {
        $cell_line_info="Missing_SAMPLE_ATTRIBUTES"
    }

$ftp_address='ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/'.substr( $sra_ID, 0, 3 ).'/'.substr( $sra_ID, 0, 6 ).'/'.$sra_ID.'/';
print OUT $sra_ID."\t".$instrument."\t".$layout."\t".$strategy."\t".$cell_line_info."\t".$organism."\t".$treated."\t".$treat_tosep."\t".$antbody."\t".$mutation1[0]."\t".$mutation2."\t".$mutations_disc."\t".$ftp_address."\t".join(';',@warning_messages)."\n";
}
close CELL_LINE_FILE;
close FACTOR_FILE;
close OUT;
exit
