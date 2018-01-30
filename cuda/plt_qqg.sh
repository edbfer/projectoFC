reset;
text_title="t = 1306"
fileout="t1306.png"
filein="dados/t1306.txt"
fileinn=""
set output
a=10.0
b=10.0
d=(2*b)/(2.0*a)
set view map
unset surface
set style data pm3d
set style function pm3d
set pm3d transparent
set noclabel
unset label
set xrange [-a:a]
set yrange [-b:b]
set xtics font "Helvetica-Bold,14"
set ytics font "Helvetica-Bold,14"
set cbrange [0:]
set format cb "%1.2e"
set isosamples 3
set size ratio d
set xlabel "x" font "Helvetica-Bold,14"
set ylabel "y" font "Helvetica-Bold,14"
set palette gray
set palette rgbformulae 30,-13,-23
set size 0.75
set terminal png
set out fileout
set title text_title  font "Helvetica-Bold,20"
set terminal pngcairo size 1920,1080 enhanced font 'Verdana,13'
set multiplot
set pointsize 1
splot filein using 1:2:3  notitle text_title with lines lc rgb "black"
print filein
set nomultiplot


