 # North boundary
 reset
 set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'
 set style data linespoints
 set grid
 set xlabel 'x (m)'
 set time
 set key left top
 set out './mach2d_output/mach2d_SIM_ID_main_u_north_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_north_boundary.dat' u 1:3 w lp t'u'
 set out './mach2d_output/mach2d_SIM_ID_main_v_north_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_north_boundary.dat' u 1:4 w lp t'v'
 set out './mach2d_output/mach2d_SIM_ID_main_p_north_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_north_boundary.dat' u 1:5 w lp t'p'
 set out './mach2d_output/mach2d_SIM_ID_main_T_north_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_north_boundary.dat' u 1:6 w lp t'T'
 
 # South boundary
 reset
 set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'
 set style data linespoints
 set grid
 set xlabel 'x (m)'
 set time
 set key left top
 set out './mach2d_output/mach2d_SIM_ID_main_u_south_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_south_boundary.dat' u 1:3 w lp t'u'
 set out './mach2d_output/mach2d_SIM_ID_main_v_south_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_south_boundary.dat' u 1:4 w lp t'v'
 set out './mach2d_output/mach2d_SIM_ID_main_p_south_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_south_boundary.dat' u 1:5 w lp t'p'
 set out './mach2d_output/mach2d_SIM_ID_main_T_south_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_south_boundary.dat' u 1:6 w lp t'T'
 
 # East boundary
 reset
 set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'
 set style data linespoints
 set grid
 set xlabel 'y (m)'
 set time
 set key left top
 set out './mach2d_output/mach2d_SIM_ID_main_u_east_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_east_boundary.dat' u 2:3 w lp t'u'
 set out './mach2d_output/mach2d_SIM_ID_main_v_east_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_east_boundary.dat' u 2:4 w lp t'v'
 set out './mach2d_output/mach2d_SIM_ID_main_p_east_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_east_boundary.dat' u 2:5 w lp t'p'
 set out './mach2d_output/mach2d_SIM_ID_main_T_east_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_east_boundary.dat' u 2:6 w lp t'T'
 
 # West boundary
 reset
 set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'
 set style data linespoints
 set grid
 set xlabel 'x (m)'
 set time
 set key left top
 set out './mach2d_output/mach2d_SIM_ID_main_u_west_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_west_boundary.dat' u 1:3 w lp t'u'
 set out './mach2d_output/mach2d_SIM_ID_main_v_west_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_west_boundary.dat' u 1:4 w lp t'v'
 set out './mach2d_output/mach2d_SIM_ID_main_p_west_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_west_boundary.dat' u 1:5 w lp t'p'
 set out './mach2d_output/mach2d_SIM_ID_main_T_west_boundary.png'
 plot './mach2d_output/mach2d_SIM_ID_main_fields_west_boundary.dat' u 1:6 w lp t'T'
