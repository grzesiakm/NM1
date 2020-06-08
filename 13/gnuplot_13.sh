set terminal post enhanced colour solid font 25  
reset

# # SKRYPT GENERUJE WYKRESY DLA TRZECH PLIKOW: c1.dat, c2.dat, c3.dat
# # W PLIKU c1.dat znajduja sie dwie kolumny dla 1. funkcji: n, blad bezwzgledny |c_1 - c_{1,a}|.
# # W PLIKU c2.dat znajduja sie dwie kolumny dla 2. funkcji: n, blad bezwzgledny |c_2 - c_{2,a}|:
#          - w 1. serii danych: wyniki dla k = 5.
#          - w 2. serii danych: wyniki dla k = 10.
# # W PLIKU c3.dat znajduja sie dwie kolumny dla 3. funkcji: n, blad bezwzgledny |c_3 - c_{3,a}|.


l_width = 5
p_size = 1.3
set style line 2 lt rgb "red" lw l_width pt 7 ps p_size
set style line 3 lt rgb "dark-spring-green" lw l_width pt 7 ps p_size
set style line 4 lt rgb "royalblue" lw l_width pt 7 ps p_size*1.2

g1(x) = x/(4.*x**2 + 1.)
g2(x, k) = x**k * exp(-x)

# 1. zadanie: wykres funkcji podcalkowej.
set samples 500
set output "c1.eps"
set xlabel "x"
set ylabel "g(x)"
set grid
plot [0:2] g1(x) ls 3 t "g_1(x)"

# 2. zadanie: dwa wykresy funkcji podcalkowej.
set output "c2_k5.eps"
plot [0:20] g2(x, 5) ls 4 t "g_2(x), k=5"

set output "c2_k10.eps"
plot [0:30] g2(x, 10) ls 4 t "g_2(x), k=10"

# 1. zadanie: wykres bledu.
set key outside top center horizontal
set output "c1_err.eps"
set xlabel "n"
set ylabel "|c_1 - c_{1,a}|"
set grid
set multiplot
set xtics out
set ytics out
plot "c1.dat" u 1:2 w lp ls 2 t "|c_1 - c_{1,a}|"
unset key
set size 0.7, 0.6
set origin 0.25, 0.3
unset xlabel
unset ylabel
set object 1 rect from graph 0, graph 0 to graph 1, graph 1  fc rgb "white" fs solid 0.15 noborder
plot [3:] "c1.dat" u 1:2 w lp ls 2 t "|c_1 - c_{1,a}|"
unset multiplot

# 2. zadanie: wykres bledu, k = 5.
set size 1.0
set key outside top center horizontal
set output "c2_err5.eps"
set xlabel "n"
set ylabel "|c_2 - c_{2,a}|"
set grid
set multiplot
unset object 1
plot "c2.dat" i 0 u 1:2 w lp ls 2 t "|c_2 - c_{2,a}|, k=5"
unset key
set size 0.7, 0.6
set origin 0.5, 0.6
unset xlabel
unset ylabel
set object 1 rect from graph 0, graph 0 to graph 1, graph 1 fc rgb "white" fs solid 0.15 noborder
plot [3:] "c2.dat" i 0 u 1:2 w lp ls 2 
unset multiplot

# 2. zadanie: wykres bledu, k = 10.
set size 1.0
set key outside top center horizontal
set output "c2_err10.eps"
set xlabel "n"
set ylabel "|c_2 - c_{2,a}|"
set grid
set multiplot
unset object 1
plot "c2.dat" i 1 u 1:2 w lp ls 2 t "|c_2 - c_{2,a}|, k=10"
unset key
set size 0.7, 0.6
set origin 0.75, 0.9
unset xlabel
unset ylabel
set xtics out
set ytics out
set object 1 rect from graph 0, graph 0 to graph 1, graph 1  fc rgb "white" fs solid 0.15 noborder
plot [6:] "c2.dat" i 1 u 1:2 w lp ls 2 t "|c_2 - c_{2,a}|"
unset multiplot

# 3. zadanie: wykres bledu.
set size 1.0
unset object 1
set key outside top center horizontal
set output "c3_err.eps"
set xlabel "n"
set ylabel "|c_3 - c_{3,a}|"
set grid
plot [:15] "c3.dat" i 0 u 1:2 w lp ls 2 t "|c_3 - c_{3,a}|"
