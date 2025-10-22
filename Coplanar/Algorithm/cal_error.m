function [error_my_R,error_my_t,error_P1P2L_R,error_P1P2L_t] = cal_error(T_est_all,T_est_all_P2P1L,T_w0_to_c0)
    R_true=T_w0_to_c0(1:3,1:3);
    t_true=T_w0_to_c0(1:3,4);
    %% 
    for i=1:4
        error1(i)=norm([acos( dot(R_true(1:3,1),T_est_all(1:3,1,i))),   ...
                        acos( dot(R_true(1:3,2),T_est_all(1:3,2,i))),   ...
                        acos( dot(R_true(1:3,3),T_est_all(1:3,3,i)))]);
        error2(i)=norm(t_true-T_est_all(1:3,4,i))/norm(t_true)*100;
    end
    index1 = find(error1==min(min(error1)));
    index2 = find(error2==min(min(error2)));
    
    error_my_R=min(error1);
    error_my_t=min(error2);
    %% 
    for i=1:4
        error1(i)=norm([acos( dot(R_true(1:3,1),T_est_all_P2P1L(1:3,1,i))),   ...
                        acos( dot(R_true(1:3,2),T_est_all_P2P1L(1:3,2,i))),   ...
                        acos( dot(R_true(1:3,3),T_est_all_P2P1L(1:3,3,i)))]);
        error2(i)=norm(t_true-T_est_all_P2P1L(1:3,4,i))/norm(t_true)*100;
    end
    index1 = find(error1==min(min(error1)));
    index2 = find(error2==min(min(error2)));
    
    error_P1P2L_R=min(error1);
    error_P1P2L_t=min(error2);
end

