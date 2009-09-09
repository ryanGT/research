

 ############
 # Title: [RigidLinkCurve]
 ############
 reset
 set terminal png small size 1000,1000;
 set output 'FIGURES\RigidLinkCurveQcv3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'RigidLinkCurve';
 set xlabel 'X2 [m]';
 set xrange [-2.50000e-001 :  2.50000e-001];
 set ylabel 'X3 [m]';
 set yrange [-2.50000e-001 :  2.50000e-001];
 plot 'FIGURES\RigidLinkCurveQCQ.mdt' using 2:3 title 'Curve' with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [RigidLinkCurve]
 ############
 reset
 set terminal png small size 1000,1000;
 set output 'FIGURES\RigidLinkCurveQcv2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'RigidLinkCurve';
 set xlabel 'X1 [m]';
 set xrange [ 7.20000e+000 :  7.70000e+000];
 set ylabel 'X3 [m]';
 set yrange [-2.50000e-001 :  2.50000e-001];
 plot 'FIGURES\RigidLinkCurveQCQ.mdt' using 1:3 title 'Curve' with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [RigidLinkCurve]
 ############
 reset
 set terminal png small size 1000,1000;
 set output 'FIGURES\RigidLinkCurveQcv1Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'RigidLinkCurve';
 set xlabel 'X1 [m]';
 set xrange [ 7.20000e+000 :  7.70000e+000];
 set ylabel 'X2 [m]';
 set yrange [-2.50000e-001 :  2.50000e-001];
 plot 'FIGURES\RigidLinkCurveQCQ.mdt' using 1:2 title 'Curve' with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size 1000, 750;
 set output 'FIGURES\PropertyBeam1QbxcQ.png';
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
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size 1000, 750;
 set output 'FIGURES\PropertyBeam1Qi11Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'AXIAL STIFFNESS [N]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QBQ.mdt' using 1:4 title 'S' with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size 1000, 750;
 set output 'FIGURES\PropertyBeam1Qi23Q.png';
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
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size 1000, 750;
 set output 'FIGURES\PropertyBeam1QbxkQ.png';
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
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size 1000, 750;
 set output 'FIGURES\PropertyBeam1Qk11Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'TORSIONAL STIFFNESS [N.m^2]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QSQ.mdt' using 1:4 title 'J' with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size 1000, 750;
 set output 'FIGURES\PropertyBeam1Qk23Q.png';
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
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size 1000, 750;
 set output 'FIGURES\PropertyBeam1QbxmQ.png';
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
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size 1000, 750;
 set output 'FIGURES\PropertyBeam1Qm00Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'PropertyBeam1';
 set xlabel 'ETA';
 set autoscale x;
 set ylabel 'MASS / SPAN [Kg/m]';
 set autoscale y;
 plot 'FIGURES\PropertyBeam1QIQ.mdt' using 1:4 title 'm00' with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [PropertyBeam1]
 ############
 reset
 set terminal png small size 1000, 750;
 set output 'FIGURES\PropertyBeam1Qm23Q.png';
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
 set output;

 ############
 # Title: [CurveBeam1]
 ############
 reset
 set terminal png small size 1000,1000;
 set output 'FIGURES\CurveBeam1Qcv3Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'CurveBeam1';
 set xlabel 'X2 [m]';
 set xrange [-3.60000e+000 :  3.60000e+000];
 set ylabel 'X3 [m]';
 set yrange [-3.60000e+000 :  3.60000e+000];
 plot 'FIGURES\CurveBeam1QCQ.mdt' using 2:3 title 'Curve' with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveBeam1]
 ############
 reset
 set terminal png small size 1000,1000;
 set output 'FIGURES\CurveBeam1Qcv2Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'CurveBeam1';
 set xlabel 'X1 [m]';
 set xrange [ 0.00000e+000 :  7.20000e+000];
 set ylabel 'X3 [m]';
 set yrange [-3.60000e+000 :  3.60000e+000];
 plot 'FIGURES\CurveBeam1QCQ.mdt' using 1:3 title 'Curve' with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveBeam1]
 ############
 reset
 set terminal png small size 1000,1000;
 set output 'FIGURES\CurveBeam1Qcv1Q.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'CurveBeam1';
 set xlabel 'X1 [m]';
 set xrange [ 0.00000e+000 :  7.20000e+000];
 set ylabel 'X2 [m]';
 set yrange [-3.60000e+000 :  3.60000e+000];
 plot 'FIGURES\CurveBeam1QCQ.mdt' using 1:2 title 'Curve' with lines linestyle 1 ;
 set nomultiplot;
 set output;

 ############
 # Title: [CurveBeam1]
 ############
 reset
 set terminal png small size 1000, 750;
 set output 'FIGURES\CurveBeam1QcrvQ.png';
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
 set output;

 ############
 # Title: [ScheduleRotationA]
 ############
 reset
 set terminal png small size 1000, 750;
 set output 'FIGURES\ScheduleRotationA.png';
 set key; set border; set grid; set nomultiplot;
 set style line  1 lt  1 lw 3.0 pt  1 ps 1.5;
 set title 'ScheduleRotationA';
 set xlabel 'TIME [sec]';
 set autoscale x;
 set ylabel 'TIME FUNCTION [sec]';
 set autoscale y;
 plot 'FIGURES\ScheduleRotationA.mdt' using 1:2 title 'History' with lines linestyle 1 ;
 set nomultiplot;
 set output;