function [T_est_all_my,time] = MV_P1P2L(delta_R,delta_t,P1_0,P2_0,L3_0,L4_0,R_w0_to_c0,t_w0_to_c0,T_w0_to_c0,noise,length)
    %
    P1_T=delta_R*R_w0_to_c0*P1_0+delta_R*t_w0_to_c0+delta_t;    P1_T=P1_T/P1_T(3);
    P1_hat=[P1_T(1)+noise/length*randn(1);P1_T(2)+noise/length*randn(1);1];
    d_1=(P1_hat)/norm(P1_hat);
    %
    P2_T=delta_R*R_w0_to_c0*P2_0+delta_R*t_w0_to_c0+delta_t;    P2_T=P2_T/P2_T(3);
    P2_hat=[P2_T(1)+noise/length*randn(1);P2_T(2)+noise/length*randn(1);1];
    d_2=(P2_hat)/norm(P2_hat);
    %
    L3_T=R_w0_to_c0*L3_0+t_w0_to_c0;   L3_T=L3_T/L3_T(3);
    L3_hat=[L3_T(1)+noise/length*randn(1);L3_T(2)+noise/length*randn(1);1];
    d_3=(L3_hat)/norm(L3_hat); 
    %
    L4_T=R_w0_to_c0*L4_0+t_w0_to_c0;   L4_T=L4_T/L4_T(3);
    L4_hat=[L4_T(1)+noise/length*randn(1);L4_T(2)+noise/length*randn(1);1];
    d_4=(L4_hat)/norm(L4_hat); 
    %% ------------------------------------------------------
    tic
    %% 
    [P1,P2,L_3,L_4,R_w1_to_w0,t_w1_to_w0] = world_translation(P1_0,P2_0,L3_0,L4_0);
    T_w1_to_w0=[R_w1_to_w0,t_w1_to_w0; 0 0 0 1];
    %% 
    [C,D_1,D_2,D_3,D_4,R_c0_to_c1,t_c0_to_c1] = camera_translation(d_1,d_2,d_3,d_4);
    T_c0_to_c1=[R_c0_to_c1,t_c0_to_c1;0 0 0 1];

    %%
    T_w1_to_c1=T_c0_to_c1*T_w0_to_c0*T_w1_to_w0;
    R_w1_to_c1=T_w1_to_c1(1:3,1:3);
    t_w1_to_c1=T_w1_to_c1(1:3,4);

    R=R_w1_to_c1;
    t=t_w1_to_c1;

    T_w1_to_c1=T_c0_to_c1*T_w0_to_c0*T_w1_to_w0;
    R_w1_to_c1=T_w1_to_c1(1:3,1:3);
    t_w1_to_c1=T_w1_to_c1(1:3,4);

    delta_T2=T_c0_to_c1*[delta_R,delta_t;[0 0 0 1]]/(T_c0_to_c1);  
    delta_R2=delta_T2(1:3,1:3);   delta_t2=delta_T2(1:3,4);
    %% ---------------------------------------------
    D_1_hat=delta_R2'*(D_1-C);   D_2_hat=delta_R2'*(D_2-C);  
    delat_t_hat=delta_R2'*(delta_t2+[0;0;1]);
    a_1_hat=D_1_hat(1); b_1_hat=D_1_hat(2);  c_1_hat=D_1_hat(3);  
    a_2_hat=D_2_hat(1); b_2_hat=D_2_hat(2);  c_2_hat=D_2_hat(3); 
    X_2=P2(1);   X_3=L_3(1);  Y_3=L_3(2); X_4=L_4(1);  Y_4=L_4(2); Z_4=L_4(3);
    delat_t_1_hat=delat_t_hat(1);  delat_t_2_hat=delat_t_hat(2);  delat_t_3_hat=delat_t_hat(3);
    % 
    x=[R(1,1),R(2,1),R(3,1),R(2,2),R(2,3),t(1),t(2),t(3)]';

    A=[0 0 0 0 0 -b_1_hat a_1_hat 0;
       0 0 0 0 0 0 -c_1_hat b_1_hat;
       -b_2_hat*X_2 a_2_hat*X_2 0 0 0 -b_2_hat a_2_hat 0;
       0 -c_2_hat*X_2 b_2_hat*X_2 0 0 0 -c_2_hat b_2_hat;
       0 X_3 0 Y_3 0 0 1 0;
       0 X_4 0 Y_4 Z_4 0 1 0];
    B=[b_1_hat*delat_t_1_hat-a_1_hat*delat_t_2_hat  ...
       c_1_hat*delat_t_2_hat-b_1_hat*delat_t_3_hat  ...
       b_2_hat*delat_t_1_hat-a_2_hat*delat_t_2_hat  ...
       c_2_hat*delat_t_2_hat-b_2_hat*delat_t_3_hat  ...
       0 0]';


    %% -----------------------------------------
 
    [R21_est,t_y_est] = solve_two_functions(A,X_3,X_4,Y_3,Y_4,Z_4,B,x);

    %% ----------------------------------------
    [R_all,t_all] = recover_T(A,X_3,X_4,Y_3,Y_4,Z_4,B,R21_est,t_y_est);
    %% 
    for i=1:4
        T_est_all_my(:,:,i)=(T_c0_to_c1)\[R_all(:,:,i),t_all(:,i);0 0 0 1]/(T_w1_to_w0);
    end
    
    time=toc;
end

