set term qt enh size 1366,768
set zeroax lt -1
set grid
set xl 'Orbits'
unset key
set sty da lin

set multiplot layout 4,2
set yl 'q_0'
p [:][-1:1]'CONTROL.dat' u ($1/6115):2
set yl 'Angle'
p '' u ($1/6115):9
set yl 'q_x'
p [:][-1:1]'' u ($1/6115):3
set yl 'w_x'
p '' u ($1/6115):6
set yl 'q_y'
p [:][-1:1]'' u ($1/6115):4
set yl 'w_y'
p '' u ($1/6115):7
set yl 'q_z'
p [:][-1:1]'' u ($1/6115):5
set yl 'w_z'
p '' u ($1/6115):8
unset multiplot

set term qt enh 2 size 1366,768
set key outside right
set multiplot layout 3,1
set yl 'Quaternion'
p [:][-1:1]'CONTROL.dat' u ($1/6115):2 t 'q_0','' u ($1/6115):3 t 'q_x','' u ($1/6115):4 t 'q_y','' u ($1/6115):5 t 'q_z'
set yl 'Angle (deg)'
p 'CONTROL.dat' u ($1/6115):9 not
set yl 'Moment (A m^{2})'
p 'EFFORT.dat' u ($1/6115):2 t 'M_x','' u ($1/6115):3 t 'M_y','' u ($1/6115):4 t 'M_z'
unset multiplot

set term qt enh 3 size 1366,768
set multiplot layout 2,1
set yl 'Quaternion'
p 'CONTROL.dat' u ($1/6115):2 t 'q_0','' u ($1/6115):3 t 'q_x','' u ($1/6115):4 t 'q_y','' u ($1/6115):5 t 'q_z'
set yl 'Angular Rate (rad s^{-1})'
p 'CONTROL.dat' u ($1/6115):6 t 'w_x','' u ($1/6115):7 t 'w_y','' u ($1/6115):8 t 'w_z'
unset multiplot

set term qt enh 4 size 1366,768
set multiplot layout 3,1
set yl 'B_{brf} (T)'
p 'MAG.dat' u ($1/6115):5 t 'B_x','' u ($1/6115):6 t 'B_y','' u ($1/6115):7 t 'B_z'
set yl 'Angle (deg)'
p 'CONTROL.dat' u ($1/6115):9 not
set yl 'Moment (A m^2)'
p 'EFFORT.dat' u ($1/6115):2 t 'M_x','' u ($1/6115):3 t 'M_y','' u ($1/6115):4 t 'M_z'
unset multiplot

