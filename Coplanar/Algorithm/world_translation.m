function [P1,P2,L_3,L_4,R_w1_to_w0,t_w1_to_w0] = world_translation(P1_0,P2_0,L3_0,L4_0)

    X2=norm(P2_0-P1_0);
    P1=[0  0 0]';    P2=[X2 0 0]';
    t_w=P1_0;
    R_x_w=(P2_0-P1_0)/norm(P2_0-P1_0);
    R_z_w=cross(R_x_w,L3_0-t_w);  R_z_w=R_z_w/norm(R_z_w);
    R_y_w=cross(R_z_w,R_x_w);R_y_w=R_y_w/norm(R_y_w);
    R_w=[R_x_w,R_y_w,R_z_w];

    
    L_3=R_w'*(L3_0-t_w);  
    L_4=R_w'*(L4_0-t_w);
    
    R_w1_to_w0=R_w;
    t_w1_to_w0=t_w;
end

