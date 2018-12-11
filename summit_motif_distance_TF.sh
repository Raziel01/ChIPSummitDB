#!/bin/sh



#max=`wc -l $1 | awk '{print $1}' `;
#echo $max;
#for ((i=1;i<=$max;i++)); do
        motif_name=`basename $1 .bed  | sed 's/_/\t/g' | awk '{print $2}'`;
#        factor_name=`awk -v i=$i -F"\t" '{OFS="\t"; if (NR==i) print $2}' $1  | sed 's/ /_\*-homerpeaks.bed \*_/g' `;
#	homerpeaks=`echo *_${factor_name}_*-homerpeaks.bed`;       
#	summit_names=`echo "${homerpeaks}" | sed 's/-homerpeaks/summit/g'`;
#	folder=`pwd`;
#	cd /molbio/projects/rSNP/newanalysis/jaspar_motifs/Peakplitter/;
#	for j in $summit_names; do 
#		b=`basename $j _summit.bed`;
#		c=`echo "$b" | sed 's/_/\t/g' | awk '{print $2"_"$3"_"$5}'`;
#		awk -v c=$c '{OFS="\t"; print $1,$2,$3,c,$4}' $j;
#	done |  grep -v "__" | sort -k1,1 -k2,2n | uniq  > /molbio/projects/rSNP/newanalysis/jaspar_motifs/analysis/statistics/summit_stat/sorted_files/sorted_summits/hs_${motif_name}_sorted_summits.bed
#	cd $folder;
	sort -k1,1 -k2,2n ${1} > analysis/statistics/summit_stat/sorted_files/motifs/hs_${motif_name}_sorted_consensus_motifset.bed;
#Legkozelebbi Summit-motif megekresese:



cd /molbio/projects/rSNP/newanalysis/jaspar_motifs/new_analysis_2017/;
echo "At `date` started:  ${motif_name} closestBed" 1>&2 ; 

	even_or_odd=`awk '{print $3-$2}' analysis/statistics/summit_stat/sorted_files/motifs/hs_${motif_name}_sorted_consensus_motifset.bed | sort | uniq`;
	isEvenNo=$( expr $even_or_odd % 2 ) 
	rm -rf analysis/statistics/summit_stat/sorted_files/closestBed/hs_${motif_name}_closest.bed;

for n in /molbio/projects/rSNP/newanalysis/jaspar_motifs/analysis/statistics_old/summit_stat/global_analysis/summit_pool/hs_all_transcriptionfactor_summits_filtered_2017_01_30_splitted/hs_all_transcriptionfactor_summits_filtered_2017_01_30_splitted*;
do
		 if [ $isEvenNo -eq 0 ]
			then
				closestBed -a $n -b analysis/statistics/summit_stat/sorted_files/motifs/hs_${motif_name}_sorted_consensus_motifset.bed | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {OFS="\t"; if ($11=="+") print $0,int($2-((($7+$8)/2))),abs(int($2-((($7+$8)/2)))); else print $0,int($2-(int(($7+$8)/2)+1)),abs(int($2-(int(($7+$8)/2)+1)))}' | awk '{if ($14<'50' && $14>'-50') print }' | awk '{OFS="\t"; print $0,$4"_"$5 }' | sort -k15,15n |  awk '!n[$16]++' | awk '{OFS="\t"; if ($11=="+") print $0; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14*-1,$15,$16}' | sed 's/(//g'  | sed 's/)//g' >> analysis/statistics/summit_stat/sorted_files/closestBed/hs_${motif_name}_closest.bed;

			else
	closestBed -a $n -b analysis/statistics/summit_stat/sorted_files/motifs/hs_${motif_name}_sorted_consensus_motifset.bed | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {OFS="\t"; print $0,int($2-(int(($7+$8)/2)+1)),abs(int($2-(int(($7+$8)/2)+1)))}' | awk '{if ($14<'50' && $14>'-50') print }' | awk '{OFS="\t"; print $0,$4"_"$5 }' | sort -k15,15n |  awk '!n[$16]++' | awk '{OFS="\t"; if ($11=="+") print $0; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14*-1,$15,$16}'  | sed 's/(//g'  | sed 's/)//g' >> analysis/statistics/summit_stat/sorted_files/closestBed/hs_${motif_name}_closest.bed;
			fi


done


	awk '{print $4}' analysis/statistics/summit_stat/sorted_files/closestBed/hs_${motif_name}_closest.bed | sed 's/__/\t/g' | awk '{print $2}' | sort | uniq > /molbio/projects/rSNP/newanalysis/jaspar_motifs/tmp/expl_${motif_name}.lst;
	z=`wc -l /molbio/projects/rSNP/newanalysis/jaspar_motifs/tmp/expl_${motif_name}.lst | awk '{print $1}'`;
	echo -e "exp_name\telement_num\tmedian\taverage\tstand.dev." > analysis/statistics/summit_stat/results/hs_${motif_name}_distancestat.tbl;
	for ((k=1;k<=$z;k++)); do
		zob=`awk -v k=$k '{ if (NR==k) print}' /molbio/projects/rSNP/newanalysis/jaspar_motifs/tmp/expl_${motif_name}.lst`;
		grep "$zob"  analysis/statistics/summit_stat/sorted_files/closestBed/hs_${motif_name}_closest.bed  | sort -k15,15n | awk '!n[$9]++' > analysis/statistics/summit_stat/experiments_summit_distance/hs_${zob}_${motif_name}.bed;
		median=`awk '{print $14}' analysis/statistics/summit_stat/experiments_summit_distance/hs_${zob}_${motif_name}.bed | Rscript /molbio/bin/czerik/statistical_analysis/SD_Mean_Median.R | grep Media | awk '{print $2}'`;
		avrg=`awk '{print $14}' analysis/statistics/summit_stat/experiments_summit_distance/hs_${zob}_${motif_name}.bed | Rscript /molbio/bin/czerik/statistical_analysis/SD_Mean_Median.R |  grep Mean | awk '{print $2}'`;
		elm=`wc -l analysis/statistics/summit_stat/experiments_summit_distance/hs_${zob}_${motif_name}.bed | awk '{print $1}'`;
		dev=`awk '{print $14}' analysis/statistics/summit_stat/experiments_summit_distance/hs_${zob}_${motif_name}.bed | Rscript /molbio/bin/czerik/statistical_analysis/SD_Mean_Median.R |  grep SD | awk '{print $2}'`;
		echo -e  "$zob\t$elm\t$median\t$avrg\t$dev" ; 
		done >> analysis/statistics/summit_stat/results/hs_${motif_name}_distancestat.tbl;


###DATABASE_table_perl_script

 awk '{print $10":"$11"__"$9}' analysis/statistics/summit_stat/sorted_files/closestBed/hs_${motif_name}_closest.bed  | sort -n | uniq  > analysis/statistics/summit_stat/database_tables/listfiles/${motif_name}_ID.lst;
perl ~/valami6.pl -ls analysis/statistics/summit_stat/database_tables/listfiles/${motif_name}_ID.lst -int1 analysis/statistics/summit_stat/sorted_files/closestBed/hs_${motif_name}_closest.bed > analysis/statistics/summit_stat/database_tables/${motif_name}_datatable.tbl;


#		cat analysis/statistics/summit_stat/experiments_summit_distance/hs_*_${motif_name}.bed | sort -k15,15n  > analysis/statistics/summit_stat/experiments_summit_distance/hs_all_${motif_name}.bed;
#OLd version:	awk -v zob=$zob '{print $0}' analysis/statistics/summit_stat/sorted_files/colosestBed/hs_${motif_name}_closest.bed  | sort -k15,15n  > analysis/statistics/summit_stat/experiments_summit_distance/hs_all_${motif_name}.bed;
#	median=`awk '{if ($11=="+") print $14; else print $14*-1}' analysis/statistics/summit_stat/experiments_summit_distance/hs_all_${motif_name}.bed | Rscript /molbio/bin/czerik/statistical_analysis/SD_Mean_Median.R |  grep Media | awk '{print $2}'`;
#	avrg=`awk '{if ($11=="+") print $14; else print $14*-1}' analysis/statistics/summit_stat/experiments_summit_distance/hs_all_${motif_name}.bed | Rscript /molbio/bin/czerik/statistical_analysis/SD_Mean_Median.R |  grep Mean | awk '{print $2}'`;
#	elm=`wc -l analysis/statistics/summit_stat/experiments_summit_distance/hs_all_${motif_name}.bed | awk '{print $1}'`;
#	dev=`awk '{if ($11=="+") print $14; else print $14*-1}' analysis/statistics/summit_stat/experiments_summit_distance/hs_all_${motif_name}.bed | Rscript /molbio/bin/czerik/statistical_analysis/SD_Mean_Median.R |  grep SD | awk '{print $2}'` ;
#	echo -e  "all\t$elm\t$median\t$avrg\t$dev" >> analysis/statistics/summit_stat/results/hs_${motif_name}_distancestat.tbl;
#done
