

 ############
 # Title: [SensorBeam1Displacements]
 ############
 set terminal postscript landscape enhanced color "Times New Roman" 24;
 set output 'E:\pythonscripts\rwkpython\TMM\beam_only\beam_onlySen.ps';
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'SensorBeam1Displacements';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'DISPLACEMENTS [m]';
 set autoscale y;
 plot 'FIGURES\SensorBeam1Displacements.mdt' using 1:2 title 'x_1' with lines linestyle 1 , \
      'FIGURES\SensorBeam1Displacements.mdt' using 1:3 title 'x_2' with lines linestyle 2 , \
      'FIGURES\SensorBeam1Displacements.mdt' using 1:4 title 'x_3' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [SensorBeam1Displacements]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'SensorBeam1Displacements';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'ROTATIONS';
 set autoscale y;
 plot 'FIGURES\SensorBeam1Displacements.mdt' using 1:5 title 'R_1' with lines linestyle 1 , \
      'FIGURES\SensorBeam1Displacements.mdt' using 1:6 title 'R_2' with lines linestyle 2 , \
      'FIGURES\SensorBeam1Displacements.mdt' using 1:7 title 'R_3' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [SensorBeam1Positions]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'SensorBeam1Positions';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'POSITION [m]';
 set autoscale y;
 plot 'FIGURES\SensorBeam1Positions.mdt' using 1:2 title 'x_1' with lines linestyle 1 , \
      'FIGURES\SensorBeam1Positions.mdt' using 1:3 title 'x_2' with lines linestyle 2 , \
      'FIGURES\SensorBeam1Positions.mdt' using 1:4 title 'x_3' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [SensorBeam1Positions]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'SensorBeam1Positions';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'ORIENTATION';
 set autoscale y;
 plot 'FIGURES\SensorBeam1Positions.mdt' using 1:5 title 'R_1' with lines linestyle 1 , \
      'FIGURES\SensorBeam1Positions.mdt' using 1:6 title 'R_2' with lines linestyle 2 , \
      'FIGURES\SensorBeam1Positions.mdt' using 1:7 title 'R_3' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [SensorBeam1Velocities]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'SensorBeam1Velocities';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'VELOCITIES [m/sec]';
 set autoscale y;
 plot 'FIGURES\SensorBeam1Velocities.mdt' using 1:2 title 'V_1' with lines linestyle 1 , \
      'FIGURES\SensorBeam1Velocities.mdt' using 1:3 title 'V_2' with lines linestyle 2 , \
      'FIGURES\SensorBeam1Velocities.mdt' using 1:4 title 'V_3' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [SensorBeam1Velocities]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'SensorBeam1Velocities';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'ANGULAR VELOCITIES [rad/sec]';
 set autoscale y;
 plot 'FIGURES\SensorBeam1Velocities.mdt' using 1:5 title 'W_1' with lines linestyle 1 , \
      'FIGURES\SensorBeam1Velocities.mdt' using 1:6 title 'W_2' with lines linestyle 2 , \
      'FIGURES\SensorBeam1Velocities.mdt' using 1:7 title 'W_3' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [SensorBeam1Forces]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'SensorBeam1Forces';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'FORCES [N]';
 set autoscale y;
 plot 'FIGURES\SensorBeam1Forces.mdt' using 1:2 title 'F_1' with lines linestyle 1 , \
      'FIGURES\SensorBeam1Forces.mdt' using 1:3 title 'F_2' with lines linestyle 2 , \
      'FIGURES\SensorBeam1Forces.mdt' using 1:4 title 'F_3' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [SensorBeam1Forces]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'SensorBeam1Forces';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'MOMENTS [N.m]';
 set autoscale y;
 plot 'FIGURES\SensorBeam1Forces.mdt' using 1:5 title 'M_1' with lines linestyle 1 , \
      'FIGURES\SensorBeam1Forces.mdt' using 1:6 title 'M_2' with lines linestyle 2 , \
      'FIGURES\SensorBeam1Forces.mdt' using 1:7 title 'M_3' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [SensorEnergy]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set title 'SensorEnergy';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'SYSTEM ENERGIES [J]';
 set autoscale y;
 plot 'FIGURES\SensorEnergy.mdt' using 1:2 title 'KeTf' with lines linestyle 1 , \
      'FIGURES\SensorEnergy.mdt' using 1:3 title 'SeTf' with lines linestyle 2 ;
 set nomultiplot;

 ############
 # Title: [SensorEnergy]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'SensorEnergy';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'TOTAL WORK [J]';
 set autoscale y;
 plot 'FIGURES\SensorEnergy.mdt' using 1:4 title 'Work' with lines linestyle 1 ;
 set nomultiplot;

 ############
 # Title: [SensorEnergy]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set title 'SensorEnergy';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'SYSTEM ENERGY JUMPS [J]';
 set autoscale y;
 plot 'FIGURES\SensorEnergy.mdt' using 1:5 title 'KeDt' with lines linestyle 1 , \
      'FIGURES\SensorEnergy.mdt' using 1:6 title 'SeDt' with lines linestyle 2 ;
 set nomultiplot;

 ############
 # Title: [SensorEnergy]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'SensorEnergy';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'TIME STEP SIZE [sec]';
 set autoscale y;
 plot 'FIGURES\SensorEnergy.mdt' using 1:7 title 'Dt' with lines linestyle 1 ;
 set nomultiplot;

 ############
 # Title: [SigPostBeam1Forces]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set nomultiplot;
 set style fill  solid 0.25 border;
 set boxwidth 0.5; set xtics 1.0;
 set title 'SigPostBeam1Forces';
 set xlabel 'NUMBER';
 set autoscale x;
 set ylabel 'MINIMUM AND MAXIMUM';
 set autoscale y;
 plot 'FIGURES\SigPostBeam1Forces.mdt' using 1:2 title 'Nb' with boxes ;
 set nomultiplot;