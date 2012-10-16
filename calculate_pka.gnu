set term gif size 1280,1024 enhanced font 'calibri,20'
set output 'titration_curve.gif'

a = 1e-5
hill = 1.
f(x) = 1 / (10 ** (hill*(a-x)) + 1)

fit f(x) 'titration_curve.dat' u 1:2 via a

set title "Titration Curve" font 'calibri,25'
set xl "pH" font 'calibri,20'
set yl "Fraction Protonated" font 'calibri,20'
#unset key

plot 'titration_curve.dat' w p pt 4 ps 2 title '', f(x) w l lw 2 lt -1 title sprintf('pK_a = %.4f',-log10(a))

