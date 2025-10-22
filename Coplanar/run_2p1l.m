
clear
clc; 
close all; clear 
Noise_Level=0;
for  monte=1:100

    [adapter,R,t,l1,l2,P1,P2,   T_point, LP1W,LP2W,T_line] = GetData(Noise_Level);
    %% 
    addpath('Algorithm')
    T_w0_to_c0=T_line;
    R_w0_to_c0=T_w0_to_c0(1:3,1:3);   t_w0_to_c0=T_w0_to_c0(1:3,4);
    delta_R=T_point(1:3,1:3)*R_w0_to_c0';
    delta_t=T_point(1:3,4)-delta_R*t_w0_to_c0;

%     tic
    [T_est_all_my,time] = MV_P1P2L(delta_R,delta_t,P1,P2,LP1W,LP2W,R_w0_to_c0,t_w0_to_c0,T_w0_to_c0,Noise_Level,1000);
    cost_time_my(monte)=time;
    [error_my_R(monte),error_my_t(monte),~,~] = ...
                      cal_error(T_est_all_my,T_est_all_my,T_w0_to_c0);

end

