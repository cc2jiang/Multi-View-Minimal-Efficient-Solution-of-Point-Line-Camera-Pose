function [C,D_1,D_2,D_3,D_4,R_c0_to_c1,t_c0_to_c1] = camera_translation(d_1,d_2,d_3,d_4)


    %
    C=[0 0 -1]';  D_3=[0 0 0 ]';

    t=C;
    R_z=d_3;
    R_y=cross(d_3,d_4);R_y=R_y/norm(R_y);
    R_x=cross(R_y,R_z);
    R=[R_x';R_y';R_z'];

    R'*(D_3-t)/norm(R'*(D_3-t));
    ddddd4=R*d_4;  a4=ddddd4(1)*1/ddddd4(3);
    D_4=[a4,0,0]';
    R'*(D_4-t)/norm(R'*(D_4-t));

    ddddd1=R*d_1;  a1=ddddd1(1)/ddddd1(3);    b1=ddddd1(2)/ddddd1(3); 
    D_1=[a1,b1,0]';
    R'*(D_1-t)/norm(R'*(D_1-t));

    ddddd2=R*d_2;  a2=ddddd2(1)/ddddd2(3);    b2=ddddd2(2)/ddddd2(3); 
    D_2=[a2,b2,0]';
    R'*(D_2-t)/norm(R'*(D_2-t));
    
    
    R_c0_to_c1=R;
    t_c0_to_c1=t;
    
end

