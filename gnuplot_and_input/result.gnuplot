set xrange[-7:7]
p \
  "spectral_function_0_0.dat"   u 1:(-$2) w l,\
  "spectral_function_2_2.dat"   u 1:(-$2) w l,\
  "spectral_function_4_4.dat"   u 1:(-$2) w l,\
  "spectral_function_0_2.dat"   u 1:(-$2) w l,\
  "spectral_function_0_4.dat"   u 1:(-$2) w l,\
  "spectral_function_2_4.dat"   u 1:(-$3) w l,\
  "spectral_function_0_0.dat_model"   u 1:(-$2) w lp,\
  "spectral_function_2_2.dat_model"   u 1:(-$2) w lp,\
  "spectral_function_4_4.dat_model"   u 1:(-$2) w lp,\
  "spectral_function_0_2.dat_model"   u 1:(-$2) w lp,\
  "spectral_function_0_4.dat_model"   u 1:(-$2) w lp,\
  "spectral_function_2_4.dat_model"   u 1:(-$3) w lp,\


pause -1

