#!/bin/bash

psffile=$1
dcdfile=$2
startframe=$3
endframe=$4

deltai=10
echo "
mol load psf $psffile 
mol addfile $dcdfile waitfor all
" > vmd_analyze_rdf.tcl
echo "Done Input"
echo " 
set sel [atomselect top \"all\"]
set n [molinfo top get numframes]
set outname \"gofr_\"
for { set i $startframe } { \$i < $endframe } { incr i $deltai} {

	\$sel frame \$i
	# \$sel update
	set gr0 [measure gofr \$sel \$sel delta 0.1 rmax 10.0 usepbc 1]
	set name [format \"%s%d.dat\" \$outname \$i]
	set outfile1 [open \$name w]

	set r [lindex \$gr0 0]
	set gr [lindex \$gr0 1]
	foreach j \$r k \$gr {
		puts \$outfile1 [format \"%.4f\t%.4f\" \$j \$k]
	}

	close \$outfile1
	# \$sel delete	 
}
exit
" >> vmd_analyze_rdf.tcl

echo " -------------------

done Input

-----------------

"


echo "
set terminal gif animate delay 20
set output '${psffile}.gif'

set xrange[0:10]
set yrange[0:3]
set style data lines
do for [i =${startframe}:${endframe}:${deltai}] { plot sprintf('gofr_%d.dat', i) using 1:2; pause 0.5 }

" > make_anim.gnu

echo "-------------------

done make_anim

-----------------
"

vmd -dispdev text -e vmd_analyze_rdf.tcl
sleep 30
gnuplot make_anim.gnu