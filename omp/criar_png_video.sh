#!/bin/bash

#Needs: gnuplot, convert and mencoder installed
#Name file with data: Density_10.dat ... Density_50.dat
#change: for r1 in {start..end..step}; do
#res: path where the data files are
#pathfig: path to move the png figures

namevid="cuda.avi"
file_ini="t"
file_end=".txt"
file_endpng=".png"
res="dados/"
pathfig="figures/"
a=10.0   #x range +- a
b=10.0   #y range +- b
nametot=""

#for r1 in {start..end..step}; do
for r1 in {1..3278..1}; do

name=$res$file_ini$r1$file_end
name0=$(basename $name .txt)
fileout=$name0".png"
filein=$name
text_title=\""t = $r1"\"

echo "reset;" > plt_qqg.sh
echo "text_title=$text_title" >> plt_qqg.sh
echo "fileout=\"$fileout\"" >> plt_qqg.sh
echo "filein=\"$filein\"" >> plt_qqg.sh
echo "fileinn=\"$filein1\"" >> plt_qqg.sh
echo "set output" >> plt_qqg.sh

echo "a=$a" >> plt_qqg.sh
echo "b=$b" >> plt_qqg.sh
echo "d=(2*b)/(2.0*a)" >> plt_qqg.sh

echo "set view map" >> plt_qqg.sh
echo "unset surface" >> plt_qqg.sh
echo "set style data pm3d" >> plt_qqg.sh
echo "set style function pm3d" >> plt_qqg.sh
echo "set pm3d transparent" >> plt_qqg.sh
echo "set noclabel" >> plt_qqg.sh
echo "unset label" >> plt_qqg.sh
echo "set xrange [-a:a]" >> plt_qqg.sh
echo "set yrange [-b:b]" >> plt_qqg.sh
echo "set xtics font \"Helvetica-Bold,14\"" >> plt_qqg.sh
echo "set ytics font \"Helvetica-Bold,14\"" >> plt_qqg.sh

echo "set cbrange [0:]" >> plt_qqg.sh
echo "set format cb \"%1.2e\"" >> plt_qqg.sh
echo "set isosamples 3" >> plt_qqg.sh
echo "set size ratio d" >> plt_qqg.sh
echo "set xlabel \"x\" font \"Helvetica-Bold,14\"" >> plt_qqg.sh
echo "set ylabel \"y\" font \"Helvetica-Bold,14\"" >> plt_qqg.sh

# FOR Gray
echo "set palette gray" >> plt_qqg.sh
# FOR Color
echo "set palette rgbformulae 30,-13,-23" >> plt_qqg.sh

echo "set size 0.75" >> plt_qqg.sh
#echo "set terminal postscript eps color enhanced" >> plt_qqg.sh
echo "set terminal png" >> plt_qqg.sh
echo "set out fileout" >> plt_qqg.sh
echo "set title text_title  font \"Helvetica-Bold,20\"" >> plt_qqg.sh
echo "set multiplot" >> plt_qqg.sh
echo "set pointsize 1"  >> plt_qqg.sh
echo "splot filein using 1:2:3  notitle text_title with lines lc rgb \"black\"" >> plt_qqg.sh
echo "print filein" >> plt_qqg.sh
echo "set nomultiplot" >> plt_qqg.sh

echo "" >> plt_qqg.sh
echo "" >> plt_qqg.sh



gnuplot plt_qqg.sh

filename0=$file_ini$r1$file_endpng
	nametot=$nametot"$pathfig$filename0,"

	echo "Trim image: $filename0"
	mv $filename0 $pathfig$filename0".png"
	convert -trim $pathfig$filename0".png" $pathfig$filename0
	rm -f $pathfig$filename0".png"
done


echo "_______________________________________________________"
echo $nametot
echo "_______________________________________________________"
echo "Creating video $namevid"
mencoder mf://$nametot -mf w=400:h=300:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o $namevid
