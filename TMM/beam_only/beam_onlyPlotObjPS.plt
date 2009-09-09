

 ############
 # Title: [PropertyBeam1]
 ############
 set terminal postscript landscape enhanced color "Times New Roman" 24;
 set output 'E:\pythonscripts\rwkpython\TMM\beam_only\beam_onlyObj.ps';
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'CENTROID LOCATION [m]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QBQ.mdt' using 1:2 title 'Xc2' with lines linestyle 1 , \
      'FIGURES\PropertyBeam1QBQ.mdt' using 1:3 title 'Xc3' with lines linestyle 2 ;
 set nomultiplot;

 ############
 # Title: [PropertyBeam1]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'AXIAL STIFFNESS [N]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QBQ.mdt' using 1:4 title 'S' with lines linestyle 1 ;
 set nomultiplot;

 ############
 # Title: [PropertyBeam1]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'BENDING STIFFNESS [N.m^2]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QBQ.mdt' using 1:5 title 'I22' with lines linestyle 1 , \
      'FIGURES\PropertyBeam1QBQ.mdt' using 1:6 title 'I33' with lines linestyle 2 , \
      'FIGURES\PropertyBeam1QBQ.mdt' using 1:7 title 'I23' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [PropertyBeam1]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'SHEAR CENTER LOCATION [m]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QSQ.mdt' using 1:2 title 'Xk2' with lines linestyle 1 , \
      'FIGURES\PropertyBeam1QSQ.mdt' using 1:3 title 'Xk3' with lines linestyle 2 ;
 set nomultiplot;

 ############
 # Title: [PropertyBeam1]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'TORSIONAL STIFFNESS [N.m^2]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QSQ.mdt' using 1:4 title 'J' with lines linestyle 1 ;
 set nomultiplot;

 ############
 # Title: [PropertyBeam1]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'SHEARING STIFFNESS [N.m^2]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QSQ.mdt' using 1:5 title 'K22' with lines linestyle 1 , \
      'FIGURES\PropertyBeam1QSQ.mdt' using 1:6 title 'K33' with lines linestyle 2 , \
      'FIGURES\PropertyBeam1QSQ.mdt' using 1:7 title 'K23' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [PropertyBeam1]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'MASS CENTER LOCATION [m]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QIQ.mdt' using 1:2 title 'Xm2' with lines linestyle 1 , \
      'FIGURES\PropertyBeam1QIQ.mdt' using 1:3 title 'Xm3' with lines linestyle 2 ;
 set nomultiplot;

 ############
 # Title: [PropertyBeam1]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'MASS / SPAN [Kg/m]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QIQ.mdt' using 1:4 title 'm00' with lines linestyle 1 ;
 set nomultiplot;

 ############
 # Title: [PropertyBeam1]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'MOMENT OF INERTIA / SPAN [Kg.m^2/m]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QIQ.mdt' using 1:5 title 'm11' with lines linestyle 1 , \
      'FIGURES\PropertyBeam1QIQ.mdt' using 1:6 title 'm22' with lines linestyle 2 , \
      'FIGURES\PropertyBeam1QIQ.mdt' using 1:7 title 'm33' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [CurveBeam1]
 ############
 set size ratio  1;
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'CurveBeam1';
 set xlabel 'X2 [m]';
 set xrange [-3.60000e+000 :  3.60000e+000];
 set ylabel 'X3 [m]';
 set yrange [-3.60000e+000 :  3.60000e+000];
 plot 'FIGURES\CurveBeam1QCQ.mdt' using 2:3 title 'Curve' with lines linestyle 1 ;
 set nomultiplot;

 ############
 # Title: [CurveBeam1]
 ############
 set size ratio  1;
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'CurveBeam1';
 set xlabel 'X1 [m]';
 set xrange [ 0.00000e+000 :  7.20000e+000];
 set ylabel 'X3 [m]';
 set yrange [-3.60000e+000 :  3.60000e+000];
 plot 'FIGURES\CurveBeam1QCQ.mdt' using 1:3 title 'Curve' with lines linestyle 1 ;
 set nomultiplot;

 ############
 # Title: [CurveBeam1]
 ############
 set size ratio  1;
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'CurveBeam1';
 set xlabel 'X1 [m]';
 set xrange [ 0.00000e+000 :  7.20000e+000];
 set ylabel 'X2 [m]';
 set yrange [-3.60000e+000 :  3.60000e+000];
 plot 'FIGURES\CurveBeam1QCQ.mdt' using 1:2 title 'Curve' with lines linestyle 1 ;
 set nomultiplot;

 ############
 # Title: [CurveBeam1]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set multiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set style line  2 lt  2 lw 3.0 pt  2 ps 1.5;
 set style line  3 lt  3 lw 3.0 pt  3 ps 1.5;
 set title 'CurveBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'CONFORMAL ROTATIONS';
 set autoscale y;
 plot 'FIGURES\CurveBeam1QTQ.mdt' using 1:2 title 'C1' with lines linestyle 1 , \
      'FIGURES\CurveBeam1QTQ.mdt' using 1:3 title 'C2' with lines linestyle 2 , \
      'FIGURES\CurveBeam1QTQ.mdt' using 1:4 title 'C3' with lines linestyle 3 ;
 set nomultiplot;

 ############
 # Title: [ScheduleRotationA]
 ############
 set size ratio 0.75;
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'ScheduleRotationA';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'TIME FUNCTION [sec]';
 set autoscale y;
 plot 'FIGURES\ScheduleRotationA.mdt' using 1:2 title 'History' with lines linestyle 1 ;
 set nomultiplot;