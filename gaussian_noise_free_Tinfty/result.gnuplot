#set term emf font "arial,22"
#set output "out.emf" 

set xtics font "arial,28"
set ytics font "arial,28"

set xlabel "{\omaga} (eV)"  font "arial,28"
set ylabel "A{\omega}(eV^{-1})"  font "arial,28"
set xrange[-7:7]



# p \
#   "Density_of_state_0_0.dat_model"  	u 1:($2) 	w lp lc rgb "black",\

p \
  "inputSpectral.dat_1_1"		u 1:(($2)) 	w lp dt 2,\
  "Density_of_state_0_0.dat"		u 1:($2) 	w lp ,\
  "Density_of_state_0_0.dat_model"  	u 1:($2) 	w lp lc rgb "black",\
  "cubic_coeff.out" u 1:2 w l lc rgb "cyan",\

pause -1
