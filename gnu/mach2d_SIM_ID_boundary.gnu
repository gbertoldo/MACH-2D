set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'
set out './mach2d_output/mach2d_SIM_ID_boundary.png'
set style data linespoints
set grid
set xlabel 'x (m)'
set ylabel 'y (m)'
set time
set size ratio -1
set key left top
set title 'Grid Boundary: SIM_ID'

plot './mach2d_output/mach2d_SIM_ID_boundary_south.dat' using 4:5 title 'south' \
    ,'./mach2d_output/mach2d_SIM_ID_boundary_north.dat' using 4:5 title 'north' \
    ,'./mach2d_output/mach2d_SIM_ID_boundary_west.dat'  using 4:5 title 'west'  \
    ,'./mach2d_output/mach2d_SIM_ID_boundary_east.dat'  using 4:5 title 'east'  

