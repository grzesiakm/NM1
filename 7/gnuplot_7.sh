set terminal post enhanced colour solid font 21  

set key spacing 2 samplen 3
set xlabel "x" 
set ylabel "f(x)"
set grid
f(x) = 1./(1.+x**2)
l_width = 5
p_size = 1.3
set samples 200
set style line 1 lt rgb "gray30" lw 6 pt 7 ps p_size
set style line 2 lt rgb "royalblue" lw l_width pt 7 ps p_size*1.2

# WYKRES DLA WEZLOW WEDLUG WZORU (2).
# ZALOZENIE: W pliku zad_1.dat znajduja sie dane w dwoch kolumnach: x, W_6(x)
set key center top
set output "zad_1.eps" 
plot [-5:5][] f(x) w l ls 1 t "f(x) = ^1/_{(1+x^2)}",\
  "zad_1.dat" w l ls 2 t "W_6(x)"

# WYKRES DLA WEZLOW ROZLOZONYCH ROWNOMIERNIE.
# ZALOZENIE: W pliku zad_2.dat znajduja sie dane w dwoch kolumnach: x, W_6(x)
set key center bottom
set output "zad_2.eps" 
plot [-5:5][] f(x) w l ls 1 t "f(x) = ^1/_{(1+x^2)}",\
  "zad_2.dat" w l ls 2 t "W_6(x)"
