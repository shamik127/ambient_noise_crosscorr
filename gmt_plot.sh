#!/bin/bash

########################################################

infile=${1}

outps="pieplots/map_receiver_array_100.ps"

rm -r $outps

gmt gmtset MAP_FRAME_TYPE fancy

gmt psbasemap -R19.7/20.1/40.5/40.8 -JM16c -B0.2/0.2WSen -K > $outps
gmt pscoast -R -J -Ia -Na -Lf-19.8/40.6/1/10+lkm -W0.5 -Dh -O -K >> $outps
gmt psscale -Ctopo.cpt -Dx8c/1.5c+w8c/0.5c+jTC+h -Bxaf -By+lsnr -O -K >> $outps
#gmt colorbar -Ct.cpt -Dx8c/1c+w12c/0.5c+jTC+h -Bxaf+l"snr" -By+lkm -O -K > $outps
#gmt makecpt -Crelief -T-8000/8000/10 -Z > grid.cpt
#gmt psimage n40_e020_1arc_v3.tif -Dx0/0+w1c+n5 -J -R -O -K >> $outps
#gmt grdimage $grd -R -J -O -K -Cgrid.cpt   >> $outps

gmt makecpt -Cjet -T0/70 -Z> topo.cpt

n=`awk '{print}' $infile | grep -v 'X2' | wc -l`
s=`seq 2 $n`


for i in $s;do
	if [[ $i -eq $n ]]
	then
		awk -F' ' -v x=$i 'NR==x {print $2,$3,$4,$5,$6,$7}' $infile | gmt psxy -R -J -Sw -Ctopo.cpt -O >> $outps
	else
		awk -F' ' -v x=$i 'NR==x {print $2,$3,$4,$5,$6,$7}' $infile | gmt psxy -R -J -Sw -Ctopo.cpt -O -K >> $outps
fi
	
done