set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'
set out './mach2d_output/mach2d_SIM_ID_residual.png'
set style data linespoints
set grid
set logscale y
set xlabel 'iteractions'
set ylabel 'norm L1 of the total residual (dimensionless)'
set format y "%2.1E"
set time
set title 'Residual of the linear systems'

plot './mach2d_output/mach2d_SIM_ID_residual.dat' using 1:2 t'LS residual' \
,'./mach2d_output/mach2d_SIM_ID_residual.dat' using 1:3  t'max(|pl|)/PF' \
,'./mach2d_output/mach2d_SIM_ID_residual.dat' using 1:12 t'mass residual' 
