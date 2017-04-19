set term qt enh size 1366,768
set grid
set zeroax lt -1
unset key
set xl 'Orbits'
set sty da lin

set multiplot layout 3,1
set yl 'B_x (T)'
p 'OUTPUT.dat' u ($1/6115):19
set yl 'B_y (T)'
p '' u ($1/6115):20
set yl 'B_z (T)'
p '' u ($1/6115):21
unset multiplot
