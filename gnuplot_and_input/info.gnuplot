#set xrange[-0:0.02]
p\
"information.out" u (log($1)):(log($2)) w lp,\
"information.out" u (log($1)):(($3)) w lp,\

#"information.out" u (1./($1)):(($3)) w lp,\


#"information.out" u (log($1)):(($2)) w lp,\
#"information.out" u (log($1)):($3) w lp,\
#"information.out" u (log($1)):($4) w lp,\
#-8.9*tanh((x-4.26)/6.6)-6.3

pause -1












#"information.out" u (1/$1):2 w lp,\
#"information.out" u (1/$1):3 w lp,\
#"information.out" u (1/$1):(($3/$1)/($2)) w lp,\
#"information.out" u (1/$1):4 w lp,\
#"information.out" u (1/$1):5 w lp


