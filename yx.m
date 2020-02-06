function yx=yx(x)
syms k A B C D;
yx=A*sin(k*x)+B*cos(k*x)+C*sinh(k*x)+D*cosh(k*x);