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

   kk(i)=vpa(fzero(@cyfc,x0(i),options));   %求k值

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
%k=1e1;    [7.2 18.1 30.4 42.6 54.7 66.9 79.1 91.3 103.5 115.6 127.8 140.0 152.2 164.3 176.5 188.7 200.9 213.0 225.2 237.4]
%k=1e2;    [7.4 18.2 30.4 42.6 54.7 66.9 79.1 91.3 103.5 115.6 127.8 140.0 152.2 164.3 176.5 188.7 200.9 213.0 225.2 237.4]
%k=1e3;    [8.9 18.3 30.4 42.6 54.8 66.9 79.1 91.3 103.5 115.6 127.8 140.0 152.2 164.3 176.5 188.7 200.9 213.0 225.2 237.4]
%k=1e4;    [12.8 19.8 30.7 42.7 54.8 67.0 79.1 91.3 103.5 115.6 127.8 140.0 152.2 164.3 176.5 188.7 200.9 213.0 225.2 237.4]
%k=1e5;    [14.9 25.5 34.2 44.0 55.4 67.3 79.3 91.4 103.5 115.7 127.9 140.0 152.2 164.4 176.5 188.7 200.9 213.1 225.2 237.4]
%k=1e6;    [15.1 27.2 39.0 50.3 61.0 71.2 81.7 92.8 104.5 116.3 128.3 140.3 152.4 164.6 176.7 188.8 201.0 213.1 225.3 237.5]
%k=1e7;    [15.2 27.3 39.5 51.6 63.7 75.7 87.6 99.4 111.1 122.6 133.9 145.0 156.2 167.5 178.9 190.6 202.4 214.2 226.2 238.2]
%k=4.5e9;    [15.219 27.397 39.574 51.750 63.927 76.103 88.279 100.456 112.632 124.807 136.983 149.158 161.333 173.508 185.683 197.857 210.031 222.205 234.378 246.551]
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
xx=0:0.001:0.258; len_w=20;

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
syms x;

E=2.1e11;
L=0.2580;
b=0.0298;
h=0.002;
In=(b*h^3)/12;
rho=7800*b*h;
m=rho;
EI=E*b*h^3/12;

global wn b1 dc;

% % % % % % % % % 撞振系统相关参数 % % % % % % % %
gk_num=2;  % 工况数目
len_w=1;   % 模态阶数 


for j=1:gk_num
    for i=1:len_w
        wn(j,i)=k(j,i)^2*sqrt(EI/m);
    end
end

% % % % % % % % % % % % % % % % % %方程中的模态质量，模态刚度参数
for j=1:gk_num
    for i=1:len_w
        phi1 = (A*sin(k(j,i)*x)+B(j,i)*cos(k(j,i)*x)-A*sinh(k(j,i)*x)-B(j,i)*cosh(k(j,i)*x))^2;
        M(j,i)=m*vpaintegral(phi1,0,L,'MaxFunctionCalls',Inf,'RelTol', 1e-32);
        phi2 = A*sin(k(j,i)*x)+B(j,i)*cos(k(j,i)*x)-A*sinh(k(j,i)*x)-B(j,i)*cosh(k(j,i)*x);
        b1(j,i)= m/M(j,i)*vpaintegral(phi2,0,L,'MaxFunctionCalls',Inf,'RelTol', 1e-32);
    end
end

%chengji=@(x,k1,B1,k2,B2) phi(x,k1,B1).*phi(x,k2,B2);
for i3=1:len_w                    %间隙-阻挡切换时的速度系数，i3代表阻挡振型
        for i4=1:len_w
            %ddy_(i3,i4)=m*(quadgk(chengji,0,L,[],[],[],[],k(2,i3),B(2,i3),k(1,i4),B(1,i4)));
            chengji = (A*sin(k(2,i3)*x)+B(2,i3)*cos(k(2,i3)*x)-A*sinh(k(2,i3)*x)-B(2,i3)*cosh(k(2,i3)*x)).*(A*sin(k(1,i4)*x)+B(1,i4)*cos(k(1,i4)*x)-A*sinh(k(1,i4)*x)-B(1,i4)*cosh(k(1,i4)*x));
            ddy_(i3,i4)=m*vpaintegral(chengji,0,L,'MaxFunctionCalls',Inf,'RelTol', 1e-32);
        end
end

for i3=1:len_w                    %间隙-阻挡切换时的速度系数，i3代表阻挡振型
        for i4=1:len_w
            %ddy_1(i4,i3)=m*(quadl(chengji,0,L,[],[],[],[],k(2,i3),B(2,i3),k(1,i4),B(1,i4)));
            chengji = (A*sin(k(2,i3)*x)+B(2,i3)*cos(k(2,i3)*x)-A*sinh(k(2,i3)*x)-B(2,i3)*cosh(k(2,i3)*x)).*(A*sin(k(1,i4)*x)+B(1,i4)*cos(k(1,i4)*x)-A*sinh(k(1,i4)*x)-B(1,i4)*cosh(k(1,i4)*x));
            ddy_1(i4,i3)=m*vpaintegral(chengji,0,L,'MaxFunctionCalls',Inf,'RelTol', 1e-32);
        end
end

for aa=1:gk_num                    %四重循环，aa,bb控制当前所处的工况振型，i,j控制当前工况下各转换时刻点到此工况下的振型
    for bb=1:len_w
        for j=1:gk_num
            for i=1:len_w
                chengji = (A*sin(k(aa,bb)*x)+B(aa,bb)*cos(k(aa,bb)*x)-A*sinh(k(aa,bb)*x)-B(aa,bb)*cosh(k(aa,bb)*x)).*(A*sin(k(j,i)*x)+B(j,i)*cos(k(j,i)*x)-A*sinh(k(j,i)*x)-B(j,i)*cosh(k(j,i)*x));
                dc(j,i,(bb-1)*gk_num+aa)=m*wn(j,i)^2* vpaintegral(chengji,0,L,'MaxFunctionCalls',Inf,'RelTol', 1e-32) /M(aa,bb);
            end
        end
    end
end

save('coefficient_deta.mat','wn','M','b1','dc','ddy_','ddy_1')
disp(['*************************************系数计算完毕**************************************']);