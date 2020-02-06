%普通悬臂梁端点基础激励振动普通数值求解，目的：为相对振型转换法做对比

clear all;clc;                      %#ok<CLALL> %
digits 21;                          %定义精度

[k,B]=modeshape;                    %求振型,输出的K B都是振型函数的参数，下面有振型函数

global w zeta_gap c1                %定义全局变量，这些全局变量在振动方程函数中要用到，具体为：ModeEqu_Gap.m中用到

w=90; zeta_gap=0.02;                %激励频率，模态阻尼取值

len_w=1;                            %模态阶数 length of mode number of w
gk_num=1;                           %gongkuang number 工况数 1为间隙内振动工况，2为接触振动工况

A=1;                                                          %以下两行定义振型函数
phi=@(x,k,B) A*sin(k*x)+B*cos(k*x)-A*sinh(k*x)-B*cosh(k*x);   %定义振型函数

L=0.2580;                           %梁的总长度
x=L;                                %定义研究点在梁上的位置
for j=1:gk_num
    for i=1:len_w
        Phi(:,i,j)=phi(x,k(j,i),B(j,i)).';                    %#ok<SAGROW> %将研究点的坐标代入振型函数，形成一个可以与模态坐标相乘的向量，这样与模态坐标相乘后就可以得到研究点的物理坐标
    end
end

load coefficient_deta.mat;          % 载入一些必要的数据，如各阶固有频率等，这些数据由Rk4sysControl_coefficient.m计算并保存到文件夹内

%                             % 计算开始
totalTime=0.5;                % 计算总时间
dt=1e-3;                      % 计算步长
timeSpan=0:dt:totalTime;      % 计算跨度

len_t=length(timeSpan);       % 计算循环的总次数
y0_=zeros(1,len_w*2);         % 数值方法需要的初值

t_=zeros(1,len_t);            % 预分配内存-时间
Eta_=zeros(1,len_w);          % 预分配内存-模态位移
y00_=zeros(len_w,1);          % 预分配内存-模态速度初值

c1=zeros(1,len_w);            % 振动方程中的常数项，每次切换方程都会改变，是相对振型转换法的精髓

eta=zeros(gk_num,len_w);      % 各个工况末端的模态置零

t0=0; t1=t0+dt;               % 首次仿真时间零点 % 首次仿真时间终点
tspan_=[t0 t1];               % 首次仿真时间历程

t_(1)=t0;                     % 仿真时刻点，对应于仿真结果的时间点 

dw(:,1)=Phi(:,:,j)*y0_(:,1);  % 端点的物理位移值

clear i x;                    % i 清除后为时间，x消除前为梁上的点位置 
j=1;                          % 工况编号 1，间隙，2，受限

for i=2:len_t
    
    t1=t0+dt;                 % 更新步长末值，步长初值在仿真后更新
    tspan_=[t0 t1];           % 更新仿真步长
        
    [tt,yy]=rk4sys(@ModeEqu_Gap,tspan_,y0_,dt);    % 仿真
    Eta_=yy(end,1);                                % 模态位移响应
      
    t_(i)=tt(end);                                % 将积分时间末值赋给输出变量
     
    dw(:,i)=Phi(:,:,1)*Eta_.'+ Phi(:,:,1)*eta(1,:).';                    % i时刻的位移响应
	eta(j,:)=eta(j,:)+Eta_;	      
	y00_=yy(end,:);                                % 记录仿真末值
    
    for i2=1:len_w
        c1(j,i2)=dc(1,:,(i2-1)*gk_num+j)*eta(1,:).';
    end
     
     y0_=[0 y00_(1,2)];     % 各阶模态位移初值置零，是相对振型转换法的精髓
                 
    t0=t_(i);                                     % 更新步长初值
    disp(['Time is ',num2str(i)]);

end

figure;
plot(t_,dw(1,:),'b');grid on;
