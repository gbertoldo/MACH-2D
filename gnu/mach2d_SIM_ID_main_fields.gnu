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

set title 'p (Pa): SIM_ID'
set out './mach2d_output/mach2d_SIM_ID_p_field.png'
splot './mach2d_output/mach2d_SIM_ID_main_fields.dat' using 1:2:3 title ''

set title 'ro (kg/m3): SIM_ID'
set out './mach2d_output/mach2d_SIM_ID_ro_field.png'
splot './mach2d_output/mach2d_SIM_ID_main_fields.dat' using 1:2:4 title ''

set title 'T (K): SIM_ID'
set out './mach2d_output/mach2d_SIM_ID_T_field.png'
splot './mach2d_output/mach2d_SIM_ID_main_fields.dat' using 1:2:5 title ''

set title 'u (m/s): SIM_ID'
set out './mach2d_output/mach2d_SIM_ID_u_field.png'
splot './mach2d_output/mach2d_SIM_ID_main_fields.dat' using 1:2:6 title ''

set title 'v (m/s): SIM_ID'
set out './mach2d_output/mach2d_SIM_ID_v_field.png'
splot './mach2d_output/mach2d_SIM_ID_main_fields.dat' using 1:2:7 title ''

set title 'M: SIM_ID'
set out './mach2d_output/mach2d_SIM_ID_M_field.png'
splot './mach2d_output/mach2d_SIM_ID_main_fields.dat' using 1:2:8 title ''

