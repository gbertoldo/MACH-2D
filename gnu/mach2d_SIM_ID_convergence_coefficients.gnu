reset
set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'
set style data pm3d
set ticslevel 0
set view map
set pm3d map
set palette
set xlabel 'x (m)'
set ylabel 'y (m)'
set time
set size ratio -1

set title 'ccu: SIM_ID'
set out './mach2d_output/mach2d_SIM_ID_ccu.png'
splot './mach2d_output/mach2d_SIM_ID_convergence_coefficients.dat' using 1:2:3 title ''

set title 'ccv: SIM_ID'
set out './mach2d_output/mach2d_SIM_ID_ccv.png'
splot './mach2d_output/mach2d_SIM_ID_convergence_coefficients.dat' using 1:2:4 title ''

set title 'ccT: SIM_ID'
set out './mach2d_output/mach2d_SIM_ID_ccT.png'
splot './mach2d_output/mach2d_SIM_ID_convergence_coefficients.dat' using 1:2:5 title ''

set title 'ccp: SIM_ID'
set out './mach2d_output/mach2d_SIM_ID_ccp.png'
splot './mach2d_output/mach2d_SIM_ID_convergence_coefficients.dat' using 1:2:6 title ''

