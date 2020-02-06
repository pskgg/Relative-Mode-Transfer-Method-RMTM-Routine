function ydot=ModeEqu_Gap(t,y)
%基础激励下 系统参数 ft=0.505cos(2*pi*fe*t)  b1为伽辽金近似后得到的参数
global zeta_gap wn b1 w c1;
% zeta_gap=0.02;
% w=120;
% 定义基础激励的参数jccc
jccc=0.505e-3*(2*pi*w)^2;
ydot=[y(2);-2*zeta_gap*wn(1,1)*y(2)-wn(1,1)^2*y(1)+b1(1,1)*jccc*sin(2*pi*w*t)-c1(1,1)];