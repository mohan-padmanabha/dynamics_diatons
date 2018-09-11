set terminal pngcairo font "arial,12"
set key top left
set output 'orientation_x.png'
set xlabel 'orientation_x'
set ylabel 'displacement'
set grid
set zeroaxis
plot 'quat_data.txt' using 8:2 title "q1" w l, \
'quat_data.txt' using 8:5 title "px" w l

set output 'orientation_y.png'
set xlabel 'orientation_y'
set ylabel 'displacement'
set grid
set zeroaxis
plot 'quat_data.txt' using 9:3 title "q2" w l, \
'quat_data.txt' using 9:6 title "py" w l



set output 'angular_vel_x.png'
set xlabel 'angular_vel_x'
set ylabel 'displacement'
set grid
set zeroaxis
plot 'vec_data.dat' using 7:1 title "ox" w l, \
'vec_data.dat' using 7:4 title "vx" w l


set output 'angular_vel_y.png'
set xlabel 'angular_vel_y'
set ylabel 'displacement'
set grid
set zeroaxis
plot 'vec_data.dat' using 8:2 title "oy" w l, \
'vec_data.dat' using 8:5 title "vy" w l


set output 'force_torque_x.png'
set xlabel 'force_torque_x'
set ylabel 'displacement'
set grid
set zeroaxis
plot 'vec_data.dat' using 7:10 title "fx" w l, \
'vec_data.dat' using 7:13 title "Tx" w l


set output 'force_torque_y.png'
set xlabel 'force_torque_y'
set ylabel 'displacement'
set grid
set zeroaxis
plot 'vec_data.dat' using 8:11 title "fy" w l, \
'vec_data.dat' using 8:14 title "Ty" w l



