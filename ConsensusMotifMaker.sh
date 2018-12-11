#!/bin/sh

#Script creating conses motif set from PWMs and ChIP-seq BED files of a scpecific transcription factor. The script includes a motif optimization step under the corresponding peak regions. The optimized motif remapped to the genome, using the identified binding sites.
cd /molbio/projects/summitdb/data_process/consensus_datatables;
if [ ! -d analysis ] ; then
	 /bin/mkdir -p analysis 
fi

if [ ! -d analysis/motif_search ] ; then
         /bin/mkdir -p analysis/motif_search
fi

if [ ! -d analysis/motif_search/homer ] ; then
         /bin/mkdir -p analysis/motif_search/homer
fi

if [ ! -d analysis/motif_search/mast ] ; then
         /bin/mkdir -p analysis/motif_search/mast
fi

if [ ! -d analysis/motif_search/fimo ] ; then
         /bin/mkdir -p analysis/motif_search/fimo
fi

if [ ! -d analysis/motif_search/combined ] ; then
         /bin/mkdir -p analysis/motif_search/combined
fi

if [ ! -d analysis/motif_search/final ] ; then
         /bin/mkdir -p analysis/motif_search/final
fi

if [ ! -d analysis/motif_search/homeropt ] ; then
         /bin/mkdir -p analysis/motif_search/homeropt
fi

if [ ! -d analysis/motif_search/homeropt/homerRes ] ; then
         /bin/mkdir -p analysis/motif_search/homeropt/homerRes
fi

if [ ! -d analysis/motif_search/homeropt/jaspar_opt_mot ] ; then
         /bin/mkdir -p analysis/motif_search/homeropt/jaspar_opt_mot
fi

max=`wc -l $1 | awk '{print $1}' `; 
echo $max;
for ((i=1;i<=$max;i++)); do
	motif_name=`awk -v i=$i -F"\t" '{OFS="\t"; if (NR==i) print $1}' $1`;
	factor_name=`awk -v i=$i -F"\t" '{OFS="\t"; if (NR==i) print $2}' $1  | sed 's/ /_\*-homerpeaks.bed \*_/g' `;
	homerpeaks=`echo *_${factor_name}_*-homerpeaks.bed`;
	##Motif optimization:
	>&2 echo "Motif optimalisation:";
	folder=`pwd`;
	cd /data/projects/summitdb/data_process/analysis/human_2018/bed; ####Loction of Homerpeaks. 
	cat $homerpeaks | sort -k1,1 -k2,2n | mergeBed | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3}' > $folder/analysis/motif_search/homeropt/hs_all_${motif_name}_peaks.bed; 
	cd `echo $folder`; 
	findMotifsGenome.pl analysis/motif_search/homeropt/hs_all_${motif_name}_peaks.bed hg19 analysis/motif_search/homeropt/homerRes/hs_all_${motif_name}_peaks -opt  /molbio/projects/rSNP/newanalysis/jaspar_motifs/homer_motif/${motif_name}_*.motif -size given;

##Creating optimal files (motif matrix and .fa files) to mast and fimo analysis:
	>&2 echo "Creating optimal files (motif matrix and .fa files) to mast and fimo analysis";
nsites=`awk '{if (NR="1") print $6}' analysis/motif_search/homeropt/homerRes/hs_all_${motif_name}_peaks/homerResults/motif1.motif | sed -e 's/^.*T://g' |  sed  's/\..*//'`;
		if [ -z "$nsites" ]; then nsites=`echo "0"`; fi
			pvalue=`awk '{if (NR="1") print $6}' analysis/motif_search/homeropt/homerRes/hs_all_${motif_name}_peaks/homerResults/motif1.motif | sed -e 's/^.*P://g' `;  
		if [ -z "$pvalue" ]; then pvalue=`echo "0"`; fi
			lenght=`wc -l analysis/motif_search/homeropt/homerRes/hs_all_${motif_name}_peaks/homerResults/motif1.motif | awk '{print $1-1}'`; # motif lenght
			tagname=`basename  /molbio/projects/rSNP/newanalysis/jaspar_motifs/homer_motif/${motif_name}_* .motif | sed 's/_/\t/g' | awk '{print $1" "$2}'`; # name and jaspar ID of the motif > needed to jaspar minimal format
		rm -rf analysis/motif_search/homeropt/jaspar_opt_mot/${motif_name}.meme; 
		echo -e "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\nMOTIF $tagname\n\nletter-probability matrix: alength= 4 w= $lenght nsites= $nsites E= $pvalue" > analysis/motif_search/homeropt/jaspar_opt_mot/${motif_name}.meme;
		awk '{OFS="\t"; if (NR>"1") print "  "$1"000",$2"000",$3"000",$4"000",""}' analysis/motif_search/homeropt/homerRes/hs_all_${motif_name}_peaks/homerResults/motif1.motif >> analysis/motif_search/homeropt/jaspar_opt_mot/${motif_name}.meme;


#bedtools getfasta -fi /home/pappcsaba/munka/hg19_ref_genome -bed analysis/motif_search/homeropt/hs_all_${motif_name}_peaks.bed -name -fo analysis/motif_search/homeropt/hs_all_${motif_name}_peaks.fa;
 
##Motif searches:
	>&2 echo "Motif search:Mast:";
#Mast_motif_search:	
	rm -rf analysis/motif_search/mast/$motif_name analysis/motif_search/mast/hs_${motif_name}_mast.bed analysis/motif_search/mast/$motif_name;
	/molbio/bin/barta/meme-old/bin/mast  analysis/motif_search/homeropt/jaspar_opt_mot/${motif_name}.meme  analysis/motif_search/homeropt/hs_all_${motif_name}_peaks.fa -oc analysis/motif_search/mast/$motif_name
	perl /molbio/projects/summitdb/data_process/consensus_datatables/new_server_scripts/scripts/MastXmlToBed_meme4.10.2_stdout.pl -in analysis/motif_search/mast/$motif_name/mast.xml  | awk '{ if ($6!="")  print  $0}'  > analysis/motif_search/mast/hs_${motif_name}_mast.bed;

# FIMO motif search
	>&2 echo "Motif search:Fimo:";
		rm -rf analysis/motif_search/fimo/hs_${motif_name}_fimo.bed analysis/motif_search/fimo/$motif_name;
		/molbio/bin/barta/meme-old/bin/fimo --max-stored-scores 10000000 --no-qvalue -o analysis/motif_search/fimo/$motif_name  analysis/motif_search/homeropt/jaspar_opt_mot/${motif_name}.meme  analysis/motif_search/homeropt/hs_all_${motif_name}_peaks.fa;
		awk -F "\t" '{OFS="\t"; gsub(/-/,"\t",$2); if (NR>"1") print $2,$3,$4,$5,$6,$7}' analysis/motif_search/fimo/$motif_name/fimo.txt |  sed 's/:/\t/g' | awk -F "\t" '{OFS="\t"; print $1,$2+$4-1,$2+$5,$1":"$2+$4"-"$3+$5,$8,$6}' > analysis/motif_search/fimo/hs_${motif_name}_fimo.bed; 
#Homer search
	>&2 echo "Motif search:Homer:";
		rm -rf analysis/motif_search/homer/hs_${motif_name}_homer.bed
	/data/programs/HOMER/bin/annotatePeaks.pl analysis/motif_search/homeropt/hs_all_${motif_name}_peaks.bed hg19 -m analysis/motif_search/homeropt/homerRes/hs_all_${motif_name}_peaks/homerResults/motif1.motif -mbed analysis/motif_search/homer/hs_${motif_name}_homer.bed -noann -nogene

	homer_min=`awk '{if (NR>"1") print }' analysis/motif_search/homer/hs_${motif_name}_homer.bed | sed 's/\./,/g' | sort -k5,5n | head -1 |  awk '{print $5+(1/1000)}'`;
	homer_max=`awk '{if (NR>"1") print }' analysis/motif_search/homer/hs_${motif_name}_homer.bed | sed 's/\./,/g' | sort -k5,5n | tail -1 |  awk '{print $5}'`;

## Creating combined motifset:
	>&2 echo "Combined motif_set";
	awk -v homer_max=$homer_max '{OFS="\t"; print $1,$2,$3,$4,homer_max,$6}'  analysis/motif_search/mast/hs_${motif_name}_mast.bed > analysis/motif_search/combined/hs_${motif_name}_combinded.bed;
	awk '{if (NR>1) print }' analysis/motif_search/homer/hs_${motif_name}_homer.bed |  intersectBed -v -a stdin -b analysis/motif_search/mast/hs_${motif_name}_mast.bed -wa -f 1  -r | intersectBed  -a stdin -b analysis/motif_search/fimo/hs_${motif_name}_fimo.bed -f 1  -r -wa >> analysis/motif_search/combined/hs_${motif_name}_combinded.bed;
	 intersectBed -v -a analysis/motif_search/fimo/hs_${motif_name}_fimo.bed -b analysis/motif_search/mast/hs_${motif_name}_mast.bed -wa -f 1  -r | intersectBed  -a stdin -b analysis/motif_search/homer/hs_${motif_name}_homer.bed -f 1  -r -wa | awk '{OFS="\t"; print $1,$2,$3,$4,"5",$6}' >> analysis/motif_search/combined/hs_${motif_name}_combinded.bed;
	 awk '{if (NR>1) print }' analysis/motif_search/homer/hs_${motif_name}_homer.bed |  intersectBed -v -a stdin -b analysis/motif_search/mast/hs_${motif_name}_mast.bed -wa -f 1  -r | intersectBed -v  -a stdin -b analysis/motif_search/fimo/hs_${motif_name}_fimo.bed -f 1  -r  >> analysis/motif_search/combined/hs_${motif_name}_combinded.bed;
	folder=`pwd`;

	
	cd `echo $folder`;

#Preparint summits and combined motif set to mtoif-summit dstance prediction
	>&2 echo "Motif summit distance calculation, creating final motif_set";
        summit_names=`echo "${homerpeaks}" | sed 's/-homerpeaks/summit/g'`;
        echo "$summit_names";
	cd /data/projects/summitdb/data_process/analysis/human_2018/bed;
        cat ${summit_names} |   awk '{if (NF==5) print}' | sort -k1,1 -k2,2n | mergeBed -d 5 | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3}' | sort -k1,1 -k2,2n  > /molbio/projects/summitdb/data_process/consensus_datatables/tmp/summits_${motif_name}.tmp;
	cd `echo $folder`;
        sort -k1,1 -k2,2n analysis/motif_search/combined/hs_${motif_name}_combinded.bed > /molbio/projects/summitdb/data_process/consensus_datatables/tmp/motif_tmp_${motif_name}.bed;
# motif-summit dstance prediction, creation consensus motif sets
	closestBed -a /molbio/projects/summitdb/data_process/consensus_datatables/tmp/summits_${motif_name}.tmp -b /molbio/projects/summitdb/data_process/consensus_datatables/tmp/motif_tmp_${motif_name}.bed | sed 's/\./,/g' | awk -F'\t' -v homer_min=$homer_min -v homer_max=$homer_max 'function abs(x){return ((x < 0.0) ? -x : x)} {OFS="\t"; if (abs((($6+$7)/2-($2+$3)/2))<30)  print $0,(($9-homer_min)/(homer_max-homer_min))*30-(abs(($6+$7)/2-($2+$3)/2)),(1+(($9-homer_min)/(homer_max-homer_min))*30)*(1/(abs((($6+$7)/2-($2+$3)/2))+1))}'  | sort -k10,10rn |  awk '!n[$4]++'  | awk '{OFS="\t"; print $5,$6,$7,$5":"$6"-"$7"_"$8,$9,$10,$11,$12}' | sort -k8,8rn |  awk '!n[$4]++' | sed 's/,/\./g' > analysis/motif_search/final/hs_${motif_name}_consensus_motifset.bed
rm -rf /molbio/projects/rSNP/newanalysis/jaspar_motifs/tmp/motif_tmp_${motif_name}.bed /molbio/projects/rSNP/newanalysis/jaspar_motifs/tmp/summits_${motif_name}.tmp; 
done
