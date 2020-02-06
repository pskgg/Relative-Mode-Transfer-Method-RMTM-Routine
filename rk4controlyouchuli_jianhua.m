%��ͨ�������˵������������ͨ��ֵ��⣬Ŀ�ģ�Ϊ�������ת�������Ա�

clear all;clc;                      %#ok<CLALL> %
digits 21;                          %���徫��

[k,B]=modeshape;                    %������,�����K B�������ͺ����Ĳ��������������ͺ���

global w zeta_gap c1                %����ȫ�ֱ�������Щȫ�ֱ������񶯷��̺�����Ҫ�õ�������Ϊ��ModeEqu_Gap.m���õ�

w=90; zeta_gap=0.02;                %����Ƶ�ʣ�ģ̬����ȡֵ

len_w=1;                            %ģ̬���� length of mode number of w
gk_num=1;                           %gongkuang number ������ 1Ϊ��϶���񶯹�����2Ϊ�Ӵ��񶯹���

A=1;                                                          %�������ж������ͺ���
phi=@(x,k,B) A*sin(k*x)+B*cos(k*x)-A*sinh(k*x)-B*cosh(k*x);   %�������ͺ���

L=0.2580;                           %�����ܳ���
x=L;                                %�����о��������ϵ�λ��
for j=1:gk_num
    for i=1:len_w
        Phi(:,i,j)=phi(x,k(j,i),B(j,i)).';                    %#ok<SAGROW> %���о��������������ͺ������γ�һ��������ģ̬������˵�������������ģ̬������˺�Ϳ��Եõ��о������������
    end
end

load coefficient_deta.mat;          % ����һЩ��Ҫ�����ݣ�����׹���Ƶ�ʵȣ���Щ������Rk4sysControl_coefficient.m���㲢���浽�ļ�����

%                             % ���㿪ʼ
totalTime=0.5;                % ������ʱ��
dt=1e-3;                      % ���㲽��
timeSpan=0:dt:totalTime;      % ������

len_t=length(timeSpan);       % ����ѭ�����ܴ���
y0_=zeros(1,len_w*2);         % ��ֵ������Ҫ�ĳ�ֵ

t_=zeros(1,len_t);            % Ԥ�����ڴ�-ʱ��
Eta_=zeros(1,len_w);          % Ԥ�����ڴ�-ģ̬λ��
y00_=zeros(len_w,1);          % Ԥ�����ڴ�-ģ̬�ٶȳ�ֵ

c1=zeros(1,len_w);            % �񶯷����еĳ����ÿ���л����̶���ı䣬���������ת�����ľ���

eta=zeros(gk_num,len_w);      % ��������ĩ�˵�ģ̬����

t0=0; t1=t0+dt;               % �״η���ʱ����� % �״η���ʱ���յ�
tspan_=[t0 t1];               % �״η���ʱ������

t_(1)=t0;                     % ����ʱ�̵㣬��Ӧ�ڷ�������ʱ��� 

dw(:,1)=Phi(:,:,j)*y0_(:,1);  % �˵������λ��ֵ

clear i x;                    % i �����Ϊʱ�䣬x����ǰΪ���ϵĵ�λ�� 
j=1;                          % ������� 1����϶��2������

for i=2:len_t
    
    t1=t0+dt;                 % ���²���ĩֵ��������ֵ�ڷ�������
    tspan_=[t0 t1];           % ���·��沽��
        
    [tt,yy]=rk4sys(@ModeEqu_Gap,tspan_,y0_,dt);    % ����
    Eta_=yy(end,1);                                % ģ̬λ����Ӧ
      
    t_(i)=tt(end);                                % ������ʱ��ĩֵ�����������
     
    dw(:,i)=Phi(:,:,1)*Eta_.'+ Phi(:,:,1)*eta(1,:).';                    % iʱ�̵�λ����Ӧ
	eta(j,:)=eta(j,:)+Eta_;	      
	y00_=yy(end,:);                                % ��¼����ĩֵ
    
    for i2=1:len_w
        c1(j,i2)=dc(1,:,(i2-1)*gk_num+j)*eta(1,:).';
    end
     
     y0_=[0 y00_(1,2)];     % ����ģ̬λ�Ƴ�ֵ���㣬���������ת�����ľ���
                 
    t0=t_(i);                                     % ���²�����ֵ
    disp(['Time is ',num2str(i)]);

end

figure;
plot(t_,dw(1,:),'b');grid on;
