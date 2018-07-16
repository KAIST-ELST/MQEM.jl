set term x11 dashed
set yrange [-4.0:2]
p\
"./reproduce_1_1.out" u 1:($3*(1))          w l  lw 3                title "1-pstar"      ,\
"./reproduce_3_3.out" u 1:($3*(1))          w l  lw 3                title "3-pstar"      ,\
"../../Sw_SOLVER.dat" u 1:($3*(1))            w lp  lw 3                title "2-Ans"      ,\
"../../Sw_SOLVER.dat" u 1:($7*(1))            w lp  lw 3                title "6-Ans"      ,\
"./reproduce_1_1.out" u 1:($2)             w l  lw 3                title "1-pstar"      ,\
"./reproduce_3_3.out" u 1:($2)             w l  lw 3                title "3-pstar"      ,\
"../Sw_SOLVER.dat" u 1:($2*(1))            w lp  lw 3                title "2-Ans"      ,\
"../Sw_SOLVER.dat" u 1:($6*(1))            w lp  lw 3                title "6-Ans"      ,\
"../Sw_SOLVER.dat" u 1:($3+$5+$7+$9+$11+$13) w l,\
-1.0/x +2.4895120644824646/x**3


pause -1
