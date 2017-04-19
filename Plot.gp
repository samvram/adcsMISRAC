set grid

set term pdfcairo enhanced
set output 'StateError.pdf'
set key outside
set zeroax lt -1

set multiplot layout 2,1 title 'State Errors'
set title 'Attitude Error'
set xlabel 'Orbits'
set ylabel 'Quaternion'
plot for [col=2:5] 'OUTPUT.dat' u ($1/6115.):col t col w l
set title 'Rate Error'
set ylabel 'Rate (rad/s)'
plot for [col=13:15] 'OUTPUT.dat' u ($1/6115.):col t col w l
unset multiplot
set output

set term pdfcairo enhanced
set output 'Requirements.pdf'

set multiplot layout 2,1 title 'Requirements'
set title 'Moment Required'
set ylabel 'Moment (Am^2)'
set xlabel 'Orbits'
plot for [col=22:24] 'OUTPUT.dat' u ($1/6115.):col t col w l
set title 'Torque Required'
set ylabel 'Torque (Nm)'
plot for [col=16:18] 'OUTPUT.dat' u ($1/6115.):col t col w l
unset multiplot
set output

set term pdfcairo enhanced
set output 'Svector.pdf'
unset key
set title 'Sliding Vector'
set xlabel 'S_x'
set ylabel 'S_y'
set zlabel 'S_z'
splot 'OUTPUT.dat' u 25:26:27 w l

set term pdfcairo
set xlabel 'x'
set ylabel 'y'
set zlabel 'z'
set output 'Orbit.pdf'
unset key
set title 'Orbit'
splot 'OUTPUT.dat' u 28:29:30 w l

set term pdfcairo
set output 'Quat.pdf'
set xlabel 'Orbits'
unset ylabel
set multiplot layout 2,2 title "Quaternion Estimation Checking"
set title 'Q0'
plot 'QUATCHECK.dat' u ($1/6115.):2 w l,'' u ($1/6115.):6 w l
set title 'Qx'
plot 'QUATCHECK.dat' u ($1/6115.):3 w l,'' u ($1/6115.):7 w l
set title 'Qy'
plot '' u ($1/6115):4 w l,'' u ($1/6115.):8 w l
set title 'Qz'
plot '' u ($1/6115):5 w l,'' u ($1/6115.):9 w l
unset multiplot

set term pdfcairo
set output 'Rate.pdf'
set xlabel 'Time (s)'
set ylabel 'Angular Rate'
set multiplot layout 2,2 title "Rate Determination Checking"
set title 'X Component'
plot 'RATECHECK.dat' u ($1/6115.):2 w l t 'Wx','' u ($1/6115.):5 w l t 'Wdx'
set title 'Y Component'
plot '' u ($1/6115.):3 w l t 'Wy','' u ($1/6115.):6 w l t 'Wdy'
set title 'Z Component'
plot '' u ($1/6115.):4 w l t 'Wz','' u ($1/6115.):7 w l t 'Wdz'
unset multiplot

