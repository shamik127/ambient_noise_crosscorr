#!/bin/bash

# ########################################################

infile=${1}

cd Results/
# awk '{print}' $infile
for d in */;do

	fl=`echo ${d} | cut -d'/' -f 1`

	km=`awk -F ','  -v y=${fl}  '{if ($1 == y) print $2} ' $infile`
	#km=`cat $infile | grep "${fl}"`
	echo $fl
	echo $km

	resultfile="${d}result.csv"
	outps="${d}pieplot.ps"
	rm -r $outps

	gmt gmtset MAP_FRAME_TYPE fancy

	gmt psbasemap -R19.7/20.1/40.5/40.8 -JM16c -B0.2/0.2WSen -K > $outps
	gmt pscoast -R -J -Ia -Na -Lf-19.8/40.6/1/10+lkm -W0.5 -Dh -O -K >> $outps
	gmt psscale -Ctopo.cpt -Dx8c/1.5c+w8c/0.5c+jTC+h -Bxaf -By+lsnr -O -K >> $outps
	gmt makecpt -Cjet -T0/${km} -Z> topo.cpt

	n=`awk '{print}' $resultfile | grep -v 'X2' | wc -l`
	s=`seq 2 $n`

	for i in $s;do
		if [[ $i -eq $n ]]
		then
			awk -F' ' -v x=$i 'NR==x {print $2,$3,$4,$5,$6,$7}' $resultfile | gmt psxy -R -J -Sw -Ctopo.cpt -O >> $outps
		else
			awk -F' ' -v x=$i 'NR==x {print $2,$3,$4,$5,$6,$7}' $resultfile | gmt psxy -R -J -Sw -Ctopo.cpt -O -K >> $outps
	fi
		
	done	
done




