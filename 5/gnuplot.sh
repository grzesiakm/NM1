set terminal png enhanced
set encoding utf8
set output  "lambda.png"
set ylabel 'λ_{k}'
set xlabel 'Numer iteracji'
plot 'lambda.dat' u 1:2 w l t 'λ_{0}' ,\
              '' u 1:3 w l t 'λ_{1}' ,\
              '' u 1:4 w l t 'λ_{2}' ,\
              '' u 1:5 w l t 'λ_{3}' ,\
              '' u 1:6 w l t 'λ_{4}' ,\
              '' u 1:7 w l t 'λ_{5}' ,\
'' u 1:8 w l t 'λ_{6}'
