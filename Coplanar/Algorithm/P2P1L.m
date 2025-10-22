function [T_est_all_P2P1L] = P2P1L(P1_0,P2_0,L3_0,L4_0,R_w0_to_c0,t_w0_to_c0,T_w0_to_c0,noise,length)
    
    %%
    %
    P1_T=R_w0_to_c0*P1_0+t_w0_to_c0;   P1_T=P1_T/P1_T(3);
    P1_hat=[P1_T(1)+noise/length*randn(1);P1_T(2)+noise/length*randn(1);1];
    d_1=(P1_hat)/norm(P1_hat); 
    %
    P2_T=R_w0_to_c0*P2_0+t_w0_to_c0;   P2_T=P2_T/P2_T(3);
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
    
    %% 坐标系转换,从W0 到W1
    [P1_T,P2_T,L_3,L_4,R_w1_to_w0,t_w1_to_w0] = world_translation(P1_0,P2_0,L3_0,L4_0);
    T_w1_to_w0=[R_w1_to_w0,t_w1_to_w0; 0 0 0 1];
    %% 坐标系转换,从C0 到C1和C1'
    [C,D_1,D_2,D_3,D_4,R_c0_to_c1,t_c0_to_c1] = camera_translation(d_1,d_2,d_3,d_4);
    T_c0_to_c1=[R_c0_to_c1,t_c0_to_c1;0 0 0 1];
    
    T_w1_to_c1=T_c0_to_c1*T_w0_to_c0*T_w1_to_w0;
    R_w1_to_c1=T_w1_to_c1(1:3,1:3);
    t_w1_to_c1=T_w1_to_c1(1:3,4);

    R=R_w1_to_c1;
    t=t_w1_to_c1;
    
    %%
    a_1=D_1(1); b_1=D_1(2);   a_2=D_2(1); b_2=D_2(2);   
    X_2=P2_T(1);   X_3=L_3(1);  Y_3=L_3(2); X_4=L_4(1);  Y_4=L_4(2); Z_4=L_4(3); 

    x=[R(1,1),R(2,1),R(3,1),R(2,2),R(2,3),t(1),t(2),t(3)]';

    A=[0 0 0 0 0 -b_1 a_1 0;
        0 0 0 0 0 0 -1 b_1;
        -b_2*X_2 a_2*X_2 0 0 0 -b_2 a_2 0;
        0 -X_2 b_2*X_2 0 0 0 -1 b_2;
        0 X_3 0 Y_3 0 0 1 0;
        0 X_4 0 Y_4 Z_4 0 1 0];
    b=[0 -b_1 0 -b_2 0 0]';
    
    A*x-b;
    %% 测试P2P1L-planar
   a1 = a_1;    b1 = b_1;
     a2 = a_2;  b2 = b_2;
     X2 = X_2;   X3 = X_3;  Y3 = Y_3;  X4 = X_4;  Y4 = Y_4; Z4 = Z_4;

    div1 = 1/(b2*a1-b1*a2);
	e1 = (b2-b1)*div1;
	e2 = (a1-a2)*div1;

	div2 = div1/Y3;
	e3 = (b1*X2*b2)*div2;
	e4 = (-X3*b2*a1 + X3*b1*a2 - b1*X2*a2)*div2;

	div3 = 1/(b1*X2*b2*div1 - Y4*e3);
    e5 = (X4+Y4*e4+b1*X2*a2*div1)*div3;
    e6 = Z4*div3;

    e7 = (e1*e5+e2);
    e8 = e1*e6;

    e9 = (e3*e5+e4);
    e10 = e3*e6;



	c1 = 1+e5*e5+e7*e7;
    c2 = 2*(e5*e6+e7*e8);
    c3 = e6*e6+e8*e8;
    c4 = 1+e9*e9;
    c5 = 2*e9*e10;
    c6 = 1+e10*e10;

	%build the final univariate quadratic equation
	f2 = (c3-c6);
	f1 = (c2-c5);
	f0 = (c1-c4);

	%solve for v
	D = f1*f1-4*f0*f2;
	if D<0
        AAAAA=1;
    end
	rD = sqrt(D);
	roots= [-(f1+rD)/(2*f2), -(f1-rD)/(2*f2)];

    for i=1:2
		v = roots(i);

		%obtain u
		u = 1/(c1+c2*v+c3*v*v);
		if u>0
            %get the rotation elements
            R21a = sqrt(u);
            R23a = R21a*v;
            R11a = e5*R21a + e6*R23a;
            R31a = e7*R21a + e8*R23a;
            R22a = e9*R21a + e10*R23a;

            %std::cout << R11a << " " << R21a << " " << R31a << " " << R22a << " " << R23a << "\n";

            T_base = X2*(a2*R21a-b2*R11a)*div1;
            T2a = b1*T_base;
            T1a = a1*T_base;
            T3a = T_base-1;
            T3b = -T_base-1;
            Ta=[T1a,T2a,T3a]';
            Tb=[-T1a,-T2a,T3b]';
            %//Deck transform: so far, everything is multiplied by -1, after that, 2 is subtracted from T3
            %//the remaining elements are calculated from two parts, the first one is multiplied by -1 and the other one is kept
            %//get the remaining elements
            div = 1/(R22a*R22a + R23a*R23a);
            R12a = (-R11a*R21a*R22a + R23a*R31a)*div;
            R13a = (-R11a*R21a*R23a - R22a*R31a)*div;
            R32a = (-R21a*R22a*R31a - R11a*R23a)*div;
            R33a = (-R21a*R23a*R31a + R11a*R22a)*div;

            r1a=[R11a,R21a,R31a]';
            r2a=[R12a,R22a,R32a]';
            r3a=[R13a,R23a,R33a]';
            
           
            
            R12b = (R11a*R21a*R22a + R23a*R31a)*div;
            R13b = (R11a*R21a*R23a - R22a*R31a)*div;
            R32b = (R21a*R22a*R31a - R11a*R23a)*div;
            R33b = (R21a*R23a*R31a + R11a*R22a)*div;


            r2b=[R12b,-R22a,R32b]';
            r3b=[R13b,-R23a,R33b]';

            %//construct the rotation matrices
 
                    
                Ra_= [r1a,r2a,r3a];
                    Rb_=[-r1a,r2b,r3b];
            %%
            R_all(:,:,2*i-1)=Ra_;
            R_all(:,:,2*i)=Rb_;
            t_all(:,2*i-1)=Ta;
            t_all(:,2*i)=Tb;
        else
            aaaaa=1;
        end
        
    end
    
%     [R_all,t_all] = recover_T_for_P1P2L(D_1,D_2,P2_T,L_3,L_4,R11_est,R21_est);
   
    %% _________________求取真正的特征到相机的姿态和位置
    for i=1:4
        T_est_all_P2P1L(:,:,i)=inv(T_c0_to_c1)*[R_all(:,:,i),t_all(:,i);0 0 0 1]*inv(T_w1_to_w0);
    end
    
end

