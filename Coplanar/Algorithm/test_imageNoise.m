clear;clc;close all
%% 世界坐标系中点特征和线特征的坐标
P1_0=0.4*[10*rand-5,10*rand-5,0]';   P2_0=0.4*[10*rand-5,10*rand-5,0]';
L3_0=0.4*[10*rand-5,10*rand-5,0]';   L4_0=0.4*[10*rand-5,10*rand-5,0]';
%% 相机坐标系，世界坐标系W0到相机坐标系C0
R_w0_to_c0=angle2dcm((180*rand-90)/180*pi,(180*rand-90)/180*pi,(180*rand-90)/180*pi,'ZYX');
R_w0_to_c0=angle2dcm(10/180*pi,20/180*pi,-20/180*pi,'ZYX');
t_w0_to_c0=[0.25*(20*rand-10),0.25*(20*rand-10),5*rand+5]';
t_w0_to_c0=[0.5,-0.5,10]';
T_w0_to_c0=[R_w0_to_c0,t_w0_to_c0; 0 0 0 1];
%% 开始蒙特卡洛仿真
data_all=[];
% 定义噪声大小
noise_all=0:0.5:2;
length=1000;
for monte_time=1:10000
    for iiii=1:size(noise_all,2)
        noise=noise_all(iiii);
        %% 一个相机开展p2p1l
        tic
        [T_est_all_P2P1L] = P2P1L(P1_0,P2_0,L3_0,L4_0,R_w0_to_c0,t_w0_to_c0,T_w0_to_c0,noise,length);
        cost_time(iiii,monte_time)=toc;
        disp(['PIP2L时间：',num2str(cost_time(iiii,monte_time))]) 
       %%
       tic
        delta_R=angle2dcm(5/180*pi,-5/180*pi,0);
        delta_t=[0,0.5,0]';
        [T_est_all_my] = MV_P1P2L(delta_R,delta_t,P1_0,P2_0,L3_0,L4_0,R_w0_to_c0,t_w0_to_c0,T_w0_to_c0,noise,length);
        cost_time_my(iiii,monte_time)=toc;
        disp(['MV-PIP2L时间：',num2str(cost_time_my(iiii,monte_time))]) 
        disp(['时间倍数：',num2str(cost_time_my(iiii,monte_time)/cost_time(iiii,monte_time))])
        %% -----------------计算姿态和位置误差
        % 经过实验验证，即时针对两点一线求出的位姿，仍然可能存在有误差为0的二义性
        [error_my_R(iiii,monte_time),error_my_t(iiii,monte_time),error_P1P2L_R(iiii,monte_time),error_P1P2L_t(iiii,monte_time)] = ...
                  cal_error(T_est_all_my,T_est_all_P2P1L,T_w0_to_c0);
    end
end
%%
meanR_P1P2L=mean(error_P1P2L_R,2);   meanT_P1P2L=mean(error_P1P2L_t,2);   mean_time_P1P2L=mean(cost_time,2);   median_time_P1P2L=median(cost_time,2);
meanR_my=mean(error_my_R,2);         meanT_my=mean(error_my_t,2);         mean_time_my=mean(cost_time_my,2); meadian_time_my=median(cost_time_my,2);
time_beishu=mean_time_my./mean_time_P1P2L;
%%
figure(1)
set(gca,'linewidth',3,'fontsize',40);
plot(noise_all,meanR_P1P2L,'k-*','linewidth',2);
hold on
plot(noise_all,meanR_my,'r-*','linewidth',2);
hold on
ylabel('Attitude error','Interpreter','latex');
% set(gca,'YLim',[70 100]);
% set(gca,'XLim',[0 200]);
xlabel('Image Noise (Pixel)','Interpreter','latex');
set(gca,'linewidth',2,'fontsize',20);
legend('P2P1L','Our Method','Location','NorthEast');
% set(h,'Interpreter','latex');
grid on

figure(2)
set(gca,'linewidth',3,'fontsize',40);
plot(noise_all,meanT_P1P2L,'k-*','linewidth',2);
hold on
plot(noise_all,meanT_my,'r-*','linewidth',2);
hold on
ylabel('Tranlation error','Interpreter','latex');
% set(gca,'YLim',[70 100]);
% set(gca,'XLim',[0 200]);
xlabel('Image Noise (Pixel)','Interpreter','latex');
set(gca,'linewidth',2,'fontsize',20);
legend('P2P1L','Our Method','Location','NorthEast');
% set(h,'Interpreter','latex');
grid on

figure(3)
set(gca,'linewidth',3,'fontsize',40);
plot(noise_all,mean_time_P1P2L,'k-*','linewidth',2);
hold on
plot(noise_all,mean_time_my,'r-*','linewidth',2);
hold on
ylabel('Time','Interpreter','latex');
% set(gca,'YLim',[70 100]);
% set(gca,'XLim',[0 200]);
xlabel('Image Noise (Pixel)','Interpreter','latex');
set(gca,'linewidth',2,'fontsize',20);
legend('P2P1L','Our Method','Location','NorthEast');
% set(h,'Interpreter','latex');
grid on







