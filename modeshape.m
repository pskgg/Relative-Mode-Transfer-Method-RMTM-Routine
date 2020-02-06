%解析法求解悬臂梁受均布荷载正弦激励端点响应
%求振型
function[k,B]=modeshape()
clear all;clc;
digits 21;

syms x k A B C D;
syms E L b H In;

bc3=E*In*d2yx(L);
bc4=E*In*d3yx(L);

bc3=expand(subs(bc3,{C,D},{-A,-B}));
bc4=expand(subs(bc4,{C,D},{-A,-B}));

d11=expand(subs(bc3,{A,B},{1,0}));
d12=expand(subs(bc3,{A,B},{0,1}));

d21=expand(subs(bc4,{A,B},{1,0}));
d22=expand(subs(bc4,{A,B},{0,1}));

cc=[d11,d12;d21,d22];
y=simplify(det(cc));

E=2.1e11;
L=0.2580;
b=0.0298;
h=0.002;
In=(b*h^3)/12;

y=subs(y);
y=subs(y,k,x);
y=eval(['@(x)',vectorize(y)]);
options=optimset('TolFun',1e-10);
x0=[7.2 18.1 30.4 42.6 54.7 66.9 79.1 91.3 103.5 115.6 127.8 140.0 152.2 164.3 176.5 188.7 200.9 213.0 225.2 237.4];

for i=1:20

   kk(i)=vpa(fzero(y,x0(i),options));   %求k值

end

%下面的程序求B值
for i=1:20
k=kk(i);
d11_=subs(d11);
d12_=subs(d12);
BB(i)=-d11_/d12_;
end

clearvars -except kk BB;

for i=1:20
        kkk(1,i) =kk(i);
		BBB(1,i) = BB(i);
end

clear kk BB;

syms x k A B C D;
syms E L b H In;
syms lambda;

bc3=E*In*d2yx(L);
bc4=E*In*d3yx(L)-lambda*yx(L);

bc3=expand(subs(bc3,{C,D},{-A,-B}));
bc4=expand(subs(bc4,{C,D},{-A,-B}));

d11=expand(subs(bc3,{A,B},{1,0}));
d12=expand(subs(bc3,{A,B},{0,1}));
d21=expand(subs(bc4,{A,B},{1,0}));
d22=expand(subs(bc4,{A,B},{0,1}));

cc=[d11,d12;d21,d22];
y=simplify(det(cc));

E=2.1e11;
L=0.2580;
b=0.0298;
h=0.002;
In=(b*h^3)/12;
lambda=1e7;

y=subs(y);
y=subs(y,k,x);
y=eval(['@(x)',vectorize(y)]);
options=optimset('TolFun',1e-10);

x0=[15.2 27.3 39.5 51.6 63.7 75.7 87.6 99.4 111.1 122.6 133.9 145.0 156.2 167.5 178.9 190.6 202.4 214.2 226.2 238.2];
%k=1e1;  [7.2 18.1 30.4 42.6 54.7 66.9 79.1 91.3 103.5 115.6 127.8 140.0 152.2 164.3 176.5 188.7 200.9 213.0 225.2 237.4]
%k=1e2;  [7.4 18.2 30.4 42.6 54.7 66.9 79.1 91.3 103.5 115.6 127.8 140.0 152.2 164.3 176.5 188.7 200.9 213.0 225.2 237.4]
%k=1e3;  [8.9 18.3 30.4 42.6 54.8 66.9 79.1 91.3 103.5 115.6 127.8 140.0 152.2 164.3 176.5 188.7 200.9 213.0 225.2 237.4]
%k=1e4;  [12.8 19.8 30.7 42.7 54.8 67.0 79.1 91.3 103.5 115.6 127.8 140.0 152.2 164.3 176.5 188.7 200.9 213.0 225.2 237.4]
%k=1e5;  [14.9 25.5 34.2 44.0 55.4 67.3 79.3 91.4 103.5 115.7 127.9 140.0 152.2 164.4 176.5 188.7 200.9 213.1 225.2 237.4]
%k=1e6;  [15.1 27.2 39.0 50.3 61.0 71.2 81.7 92.8 104.5 116.3 128.3 140.3 152.4 164.6 176.7 188.8 201.0 213.1 225.3 237.5]
%k=1e7;  [15.2 27.3 39.5 51.6 63.7 75.7 87.6 99.4 111.1 122.6 133.9 145.0 156.2 167.5 178.9 190.6 202.4 214.2 226.2 238.2]

for i=1:20

   kk(i)=vpa(fzero(y,x0(i),options));   %求k值

end

%下面的程序求B值
for i=1:20
k=kk(i);
d11_=subs(d11);
d12_=subs(d12);
BB(i)=-d11_/d12_;
end

clearvars -except kkk BBB kk BB;
%%  (1,*)表示悬臂梁参数，(2,*)表示受限制梁参数
for i=1:20
       	k(1,i) = kkk(i);
		B(1,i) = BBB(i);
        k(2,i) = kk(i);
		B(2,i) = BB(i);       
end

clear kk BB kkk BBB;

% % % % % % % % % % % % % % % % % 振型函数 % % % % % % % % % % % % % % % % % 
A=1;
phi=@(x,k,B) A*sin(k*x)+B*cos(k*x)-A*sinh(k*x)-B*cosh(k*x);

% % % % % % % % % % % % % % % % % %画振型图验证开始
xx=0:0.001:0.258; len_w=10;

for i=1:len_w
       yy= phi(xx,k(1,i),B(1,i));
	   plot(xx,yy);hold on; 
end

figure;
for i=1:len_w
	   yy= phi(xx,k(2,i),B(2,i));
	   plot(xx,yy,'r');hold on;
end
% % % % % % % % % % % % % % % % %画振型图验证结束