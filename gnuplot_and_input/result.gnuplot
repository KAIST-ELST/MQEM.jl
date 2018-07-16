set xrange[-7:7]
p \
  "Density_of_state_0_0.dat"   u 1:(-$2) w l,\
  "Density_of_state_2_2.dat"   u 1:(-$2) w l,\
  "Density_of_state_4_4.dat"   u 1:(-$2) w l,\
  "Density_of_state_0_2.dat"   u 1:(-$2) w l,\
  "Density_of_state_0_4.dat"   u 1:(-$2) w l,\
  "Density_of_state_2_4.dat"   u 1:(-$3) w l,\
  "Density_of_state_0_0.dat_model"   u 1:(-$2) w lp,\
  "Density_of_state_2_2.dat_model"   u 1:(-$2) w lp,\
  "Density_of_state_4_4.dat_model"   u 1:(-$2) w lp,\
  "Density_of_state_0_2.dat_model"   u 1:(-$2) w lp,\
  "Density_of_state_0_4.dat_model"   u 1:(-$2) w lp,\
  "Density_of_state_2_4.dat_model"   u 1:(-$3) w lp,\


pause -1

