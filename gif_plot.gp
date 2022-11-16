#! /opt/homebrew/bin/gnuplot

cd 'data'

set terminal gif animate delay 3
set output 'gif/rhoPlot.gif'
set hidden3d
set xyplane at 0
set xlabel "x"
set ylabel "y"
set zlabel "rho(x,y)"
set zrange [0.0:2.0]
set cbrange [0.0:2.0]
stats 'rhoResults.dat' nooutput
do for [i=1:int(STATS_blocks)-1] {
    splot 'rhoResults.dat' index (i-1) using 3:2:4 with lines palette notitle
}

set terminal gif animate delay 3
set output 'gif/momentumX_Plot.gif'
set hidden3d
set xyplane at 0
set xlabel "x"
set ylabel "y"
set zlabel "momentumX(x,y)"
set zrange [0.0:2.0]
set cbrange [0.0:2.0]
stats 'momentumX_Results.dat' nooutput
do for [i=1:int(STATS_blocks)-1] {
    splot 'momentumX_Results.dat' index (i-1) using 2:3:4 with lines palette notitle
}

set terminal gif animate delay 3
set output 'gif/momentumY_Plot.gif'
set hidden3d
set xyplane at 0
set xlabel "x"
set ylabel "y"
set zlabel "momentumY(x,y)"
set zrange [0.0:2.0]
set cbrange [0.0:2.0]
stats 'momentumY_Results.dat' nooutput
do for [i=1:int(STATS_blocks)-1] {
    splot 'momentumY_Results.dat' index (i-1) using 2:3:4 with lines palette notitle
}

set terminal gif animate delay 3
set output 'gif/energyPlot.gif'
set hidden3d
set xyplane at 0
set xlabel "x"
set ylabel "y"
set zlabel "energy(x,y)"
set zrange [0.0:8.0]
set cbrange [0.0:8.0]
stats 'energyResults.dat' nooutput
do for [i=1:int(STATS_blocks)-1] {
    splot 'energyResults.dat' index (i-1) using 2:3:4 with lines palette notitle
}
