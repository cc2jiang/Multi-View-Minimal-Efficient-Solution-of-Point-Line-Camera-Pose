function [adapter,R,C,l1,l2,P1,P2,   T_point,     LP1W,LP2W,T_line]=GetData(Noise_Level)


if nargin == 0
    Plot = 1;
    Noise_Level = 0;
end

K = [1000, 0.000000e+00, 6.900000e+02;
     0.000000e+00, 1000, 2.331966e+02;
     0.000000e+00, 0.000000e+00, 1.000000e+00];

dangle = 2;

 
% get rand data
C = [randn()*5-10; randn()*5; randn()*5];
a1 = randn()*pi/9; a2 = randn()*pi/9; a3 = randn()*pi/9;
R = angle2dcm(a1,a2,a3);

P1  = [randn()*2.5-15;randn()*5+15;randn()*5];
P2  = [randn()*2.5-10;randn()*5-15;randn()*5];
LP1W = [randn()*2.5;randn()*2.5-7.5;randn()*2.5];
LP2W = [randn()*2.5;randn()*2.5+7.5;randn()*2.5];
L1  = [LP1W,LP2W];

C1C = [0;5;0]+randn(3,1)*2;
C1W = R*C1C + C;
D1W = P1-C1W; D1W = D1W/norm(D1W);
D1C = R'*D1W; D1C = D1C/norm(D1C);
Rcamera1 = [1,0,0;0,0,-1;0,1,0]*angle2dcm(randn()*dangle*pi/180,randn()*dangle*pi/180,randn()*dangle*pi/180);

%%
T_point=[Rcamera1*R',-Rcamera1*R'*C1W;[0 0 0 1]];

if Noise_Level ~= 0
    
    u   = K*Rcamera1*D1C;
    u   = u/u(3);
    D1C = K\(u + [randn(2,1)*Noise_Level; 0]);
    D1C = Rcamera1'*(D1C/norm(D1C));
end

C2C = C1C;
C2W = R*C2C + C;
D2W = P2-C2W; D2W = D2W/norm(D2W);
D2C = R'*D2W; D2C = D2C/norm(D2C);
Rcamera2=Rcamera1;
if Noise_Level ~= 0
    u   = K*Rcamera2*D2C;
    u   = u/u(3);
    D2C = K\(u + [randn(2,1)*Noise_Level; 0]);
    D2C = Rcamera2'*(D2C/norm(D2C));
end

LP1C = R'*LP1W -R'*C;
LP2C = R'*LP2W -R'*C;
LP1C = LP1C/norm(LP1C);
LP2C = LP2C/norm(LP2C);
Rcamera0 = [1,0,0;0,1,0;0,0,1]*angle2dcm(randn()*dangle*pi/180,randn()*dangle*pi/180,randn()*dangle*pi/180);

%%
T_line=[Rcamera0*R',-R'*C;[0 0 0 1]];

if Noise_Level ~= 0
    u   = K*Rcamera0*LP2C;
    u   = u/u(3);
%     LP2C = K\(u + [randn(2,1)*Noise_Level/2; 0]);
    LP2C = K\(u + [randn(2,1)*Noise_Level; 0]);
    LP2C = Rcamera0'*LP2C/norm(LP2C);
    
    u   = K*Rcamera0*LP1C;
    u   = u/u(3);
%     LP1C = K\(u + [randn(2,1)*Noise_Level/2; 0]);
    LP1C = K\(u + [randn(2,1)*Noise_Level; 0]);
    LP1C = Rcamera0'*(LP1C/norm(LP1C));
end
N1C = cross(LP1C,LP2C);
N1C = N1C/norm(N1C);

% output data
l1 = norm(P1-C1W);
l2 = norm(P2-C2W);

% set the adapter
% camera info
adapter.camera.c1 = C1C;
adapter.camera.c2 = C2C;
adapter.camera.d1 = D1C;
adapter.camera.d2 = D2C;
adapter.camera.l1.normal = N1C;
% world info
adapter.world.p1 = P1;
adapter.world.p2 = P2;
adapter.world.l1.p1 = LP1W;
adapter.world.l1.p2 = LP2W;


end
