#!/bin/bash

file_ini="t"
file_endpng=".png"
namevid=gpe.avi
pathfig="figures/"
nametot=""
for r1 in {1..5000..1}; do
filename0=$file_ini$r1$file_endpng
nametot=$nametot"$pathfig$filename0,"
done

echo "_______________________________________________________"
echo $nametot
echo "_______________________________________________________"
echo "Creating video $namevid"
mencoder mf://$nametot -mf w=1920:h=1080:fps=60:type=png -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o $namevid
