#!/bin/bash

sd=1;
for ((r=1;r<=1000;r++)); do
	for m in GTR+FO+G8; do
		ncount=-1;
		for n in 100 200 300 400 500 750 1000 1500 2000; do
			if [[ ! -e ${n}bp/ ]]; then
				mkdir ${n}bp;
			fi
			cd ${n}bp;
			ncount=$((ncount+1));
			if [[ $ncount -lt 4 ]]; then
				thrds=2
			elif [[ $ncount -gt 4 && $ncount -lt 7 ]]; then
				thrds=3
			else
				thrds=4
			fi
			for g in 0.00 0.10 0.20 0.30 0.40; do
				if [[ ! -e gp${g}/ ]]; then
					mkdir gp${g};
				fi
				cd gp${g};
				for t in 30 80 130 180; do
					if [[ ! -e ${t}tips/ ]]; then
						mkdir ${t}tips;
					fi
					cd ${t}tips;
					if [[ ! -e subdir${sd}/ ]]; then
						mkdir subdir${sd};
					fi
					cd subdir${sd};
					until [[ $(ls | wc -l) -lt 1000 ]]; do
						cd ../
						sd=$((sd+1));
						if [[ ! -e subdir${sd}/ ]]; then
							mkdir subdir${sd};
						fi
						cd subdir${sd};
					done
					if [[ ! -e deltacalc.R ]]; then
						cp ../../../../*.pl .;
						cp ../../../../deltacalc.R .;
					fi
					s1=$SECONDS;
					d2=0;
					d3=0;
					if [[ ! $(grep ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}: ../../../../stopwatch.txt) ]]; then
						# verifica comprimento da árvore template
						l=$(perl treelength.pl ../../../../true.tree_${m/+FO+G8/}_${t}tips.tre);
						if [[ $(perl -e 'if($ARGV[0]>0){print 1}else{print 0}' ${g}) -eq 1 ]]; then
							# modelo preditor do indel rate via tamanho da arvore = (0.11545-(0.01669*$l)+${g})/9.69185
							# modelo preditor de indel rate via numero de tips = (0.0720572-(0.0006168*$t)+${g})/9.8042663
							i=$(perl -e 'my $i=(0.11545-(0.01669*$ARGV[0])+$ARGV[1])/9.69185; if($i>0){printf("%.4f",$i)}else{print 0.001};' ${l} ${g});
							echo "Estimated indel rate: $i";
						else
							i=0;
						fi
						dif=0;
						# formatacao do arquivo controle do indelible
						perl iqtree2indelible.pl ../../../../FoleySpringerTeelingPRSB2016_${m}_${t}taxa.iqtree ../../../../controlModel.txt ${n} INDR > control_${r}_temp.txt;
						# repete simulacao do alinhamento ate que tenha determinada proporcao de gaps \/
						s2=$SECONDS;
						mkdir sim_rep${r};
						until [[ -s simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.noGapOnly.${n}bp.phy && $(perl -e 'if($ARGV[0]<0.01){print 1}else{print 0};' ${dif}) -eq 1 ]]; do
							if [[ -s simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.noGapOnly.${n}bp.phy ]]; then
								if [[ $(perl -e 'if($ARGV[0]>$ARGV[1]){print 1}else{print 0};' ${g} ${GP}) -eq 1 ]]; then
									i=$(perl -e 'my $i=$ARGV[0]+($ARGV[1]/10); print $i;' ${i} ${dif});
								else
									i=$(perl -e 'my $i=$ARGV[0]-($ARGV[1]/10); if($i>0){print $i}else{print 0.001};' ${i} ${dif});
								fi
							fi
							sed "s/GAPPROP/$g/;s/REP/$r/g;s/INDR/$i/g" control_${r}_temp.txt > sim_rep${r}/control.txt;
							cd sim_rep${r};
							indelible;
							mv simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.phy ..;
							cd ..;
							perl noGapOnlySites.pl simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.phy > simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.noGapOnly.phy;
							perl randomSiteSampler.pl $n simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.noGapOnly.phy;
							perl fasta2phylip-vice-versa.pl simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.noGapOnly.${n}bp.phy;
							GP=$(perl gapcounter.pl simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.noGapOnly.${n}bp.fasta);
							dif=$(perl -e 'if($ARGV[0]>=$ARGV[1]){print $ARGV[0]-$ARGV[1]}else{print $ARGV[1]-$ARGV[0]};' ${g} ${GP});
							echo "Indel rate: ${i}";
							echo "|${g}-${GP}| = ${dif}";
						done
						mv simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.noGapOnly.${n}bp.phy simdata_${m/+FO+G8/}_gp${g}_${n}bp_${t}seqs_rep${r}.phy;
						rm simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.phy;
						rm simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.noGapOnly.phy;
						rm simdata_${m/+FO+G8/}_gp${g}_${n}bp_rep${r}_TRUE.noGapOnly.${n}bp.fasta;
						rm control_${r}_temp.txt;
						rm -r sim_rep${r};
						d2=$((SECONDS-s2));
						# estima parâmetros e otimiza coprimentos de ramo de cada árvore via ML com base no alinhamento simulado \/
						s3=$SECONDS;
						iqtree-omp -s simdata_${m/+FO+G8/}_gp${g}_${n}bp_${t}seqs_rep${r}.phy -te ../../../../true.tree_${m/+FO+G8/}_${t}tips.tre -pre ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r} -m ${m} -nt ${thrds} -redo -quiet;
						perl noShortLeaf.pl ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.treefile 0;
						tree="${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.treefile.noShortLeaf";
						params=$(grep -o 'parameters):.*' ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.iqtree | cut -d: -f2);
						tips=$(grep -o '[a-z]:' ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.treefile.noShortLeaf | wc -l);
						lkhd=$(grep 'BEST SCORE' ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.log | cut -d: -f2);
						gamma=$(grep 'Gamma shape alpha:' ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.log | cut -d: -f2);
						treel=$(perl treelength.pl ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.treefile.noShortLeaf);
						rm ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.ckp.gz;
						if [[ -e ${m/+FO+G8/}_gp${g}_${n}bp_${t}seqs_rep${r}.uniqueseq.phy ]]; then
							ugp=$(perl gapcounter.pl ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.uniqueseq.phy);
							rm ${m/+FO+G8/}_gp${g}_${n}bp_${t}seqs_rep${r}.uniqueseq.phy;
						else
							ugp=$(perl gapcounter.pl simdata_${m/+FO+G8/}_gp${g}_${n}bp_${t}seqs_rep${r}.phy);
						fi
						echo "${tree},${params},${tips},${lkhd},${gamma},${n},${ugp},${treel},NA,NA,NA,NA,NA,NA,NA,NA" > ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.compdatatemp
						for d in 0.05 0.10 0.15 0.20 0.25 0.30; do
							for s in 1 2 3; do
								for v in 1 2 3; do
									iqtree-omp -s simdata_${m/+FO+G8/}_gp${g}_${n}bp_${t}seqs_rep${r}.phy -te ../../../../reag.tree_${m/+FO+G8/}_${t}tips_dist${d}_set${s}_t${v}.tre -pre ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}_dist${d}_set${s}_t${v} -m ${m} -nt ${thrds} -redo -quiet;
									perl noShortLeaf.pl ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}_dist${d}_set${s}_t${v}.treefile 0;
									tree="${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}_dist${d}_set${s}_t${v}.treefile.noShortLeaf";
									lkhd=$(grep 'BEST SCORE' ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}_dist${d}_set${s}_t${v}.log | cut -d: -f2);
									gamma=$(grep 'Gamma shape alpha:' ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}_dist${d}_set${s}_t${v}.log | cut -d: -f2);
									treel=$(perl treelength.pl ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}_dist${d}_set${s}_t${v}.treefile.noShortLeaf);
									echo "${tree},${params},${tips},${lkhd},${gamma},${n},${ugp},${treel},NA,NA,NA,NA,NA,NA,NA,NA" >> ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.compdatatemp;
									rm ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}_dist${d}_set${s}_t${v}.ckp.gz;
									if [[ -e ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}_dist${d}_set${s}_t${v}.uniqueseq.phy ]]; then
										rm ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}_dist${d}_set${s}_t${v}.uniqueseq.phy;
									fi
								done
							done
						done					
						Rscript --default-packages=phangorn,utils,methods deltacalc.R < ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.compdatatemp > ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.compdata;
						cat ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.compdata >> ../../../../frame.txt;
						rm ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.compdatatemp;
						rm ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}.compdata;
						rm simdata_${m/+FO+G8/}_gp${g}_${n}bp_${t}seqs_rep${r}.phy;
						rm ${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}*.treefile.noShortLeaf;
						d3=$((SECONDS-s3));
						dt=$((SECONDS-s1));
						echo "${m/+FO+G8/}_gp${g}_${n}bp_${t}tips_rep${r}: ${d2}s ${d3}s ${dt}s" >> ../../../../stopwatch.txt;
					fi
					cd ..;
					cd ..;
				done
				cd ..;
			done
			cd ..;
		done
	done
done
