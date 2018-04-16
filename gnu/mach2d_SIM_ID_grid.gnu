reset
set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'
set out './mach2d_output/mach2d_SIM_ID_grid.png'
set style data lines
set grid
set xlabel 'x (m)'
set ylabel 'y (m)'
set time
set size ratio -1
set title 'Grid: SIM_ID'
plot './mach2d_output/mach2d_SIM_ID_grid.dat' using 1:2 title ''



