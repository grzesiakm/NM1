set terminal post enhanced colour solid font 21 

set key spacing 2 samplen 3
set samples 200
set xlabel "x" 
set ylabel "f(x)"
set grid
f1(x) = 1./(1.+x**2)
f2(x) = cos(2.*x)
l_width = 5
p_size = 1.3

set style line 1 lt rgb "gray30" lw 6 pt 7 ps p_size
set style line 2 lt rgb "royalblue" lw l_width pt 7 ps p_size*1.2

# 4 WYKRESY DLA FUNKCJI f1.
# ZALOZENIE: W pliku f1.dat znajduja sie dane w dwoch kolumnach (x oraz s(x)). 
# Pierwsza seria danych: dla n=5, nastepnie dwie puste linie,
# druga seria danych (dla n=6), dwie puste linie, itd.

i=0; n=5
set output "f1_n".n.".eps" 
plot [-5:5][] f1(x) w l ls 1 t "f_1(x) = ^1/_{(1+x^2)}",\
  "f1.dat" i i u 1:2 w l ls 2 t "s(x), n = ".n.""

i=1; n=6
set output "f1_n".n.".eps" 
plot [-5:5][] f1(x) w l ls 1 t "f_1(x) = ^1/_{(1+x^2)}",\
  "f1.dat" i i u 1:2 w l ls 2 t "s(x), n = ".n.""

i=2; n=10
set output "f1_n".n.".eps" 
plot [-5:5][] f1(x) w l ls 1 t "f_1(x) = ^1/_{(1+x^2)}",\
  "f1.dat" i i u 1:2 w l ls 2 t "s(x), n = ".n.""

i=3; n=20
set output "f1_n".n.".eps" 
plot [-5:5][] f1(x) w l ls 1 t "f_1(x) = ^1/_{(1+x^2)}",\
  "f1.dat" i i u 1:2 w l ls 2 t "s(x), n = ".n.""


# 3 WYKRESY DLA FUNKCJI f2.
# ZALOZENIE: W pliku f2.dat znajduja sie dane w dwoch kolumnach (x oraz s(x)). 
# Pierwsza seria danych: dla n=6, nastepnie dwie puste linie,
# druga seria danych (dla n=7), dwie puste linie, itd.

set ytics 0.5
set key right top horizontal outside

i=0; n=6
set output "f2_n".n.".eps" 
plot [-5:5][] f2(x) w l ls 1 t "f_2(x) = cos(2x)",\
  "f2.dat" i i u 1:2 w l ls 2 t "s(x), n = ".n.""

i=1; n=7
set output "f2_n".n.".eps" 
plot [-5:5][] f2(x) w l ls 1 t "f_2(x) = cos(2x)",\
  "f2.dat" i i u 1:2 w l ls 2 t "s(x), n = ".n.""

i=2; n=14
set output "f2_n".n.".eps" 
plot [-5:5][] f2(x) w l ls 1 t "f_2(x) = cos(2x)",\
  "f2.dat" i i u 1:2 w l ls 2 t "s(x), n = ".n.""

