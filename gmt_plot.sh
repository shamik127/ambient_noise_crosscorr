#!/bin/bash

########################################################

infile=${1}
grd=./data.grd 

outps="map_receiver_array.ps"

rm -r $outps

gmt gmtset MAP_FRAME_TYPE fancy
gmt psbasemap -R19.7/20.1/40.5/40.8 -JM16c -B0.2/0.2WSen -K > $outps

gmt makecpt -C#721b65,#b80d57,#f8615a,#ffd868 -T0/1 -Z > topo.cpt
#gmt pscoast -R -J -Ia -Na -Lf-19.8/40.6/1/10+lkm -W0.5 -Dh -O -K >> $outps
#gmt makecpt -Crelief -T-8000/8000/10 -Z > grid.cpt
gmt psimage n40_e020_1arc_v3.tif -Dx0/0+w1c+n5 -J -R -O -K >> $outps

#gmt grdimage $grd -R -J -O -K -Cgrid.cpt   >> $outps


n=`awk '{print}' result.csv | grep -v 'X2' | wc -l`
s=`seq 2 $n`


for i in $s;do
	if [[ $i -eq $n ]]
	then
		awk -F',' -v x=$i 'NR==x {print $2,$3,$4,$5,$6,$7}' $infile | gmt psxy -R -J -Sw -Ctopo.cpt -O >> $outps
	else
		awk -F',' -v x=$i 'NR==x {print $2,$3,$4,$5,$6,$7}' $infile | gmt psxy -R -J -Sw -Ctopo.cpt -O -K >> $outps
fi
	
done
#gv $outps &