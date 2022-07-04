 function [M01_,M1_1,M12,M23] = TR4_positionMat(Q,L)
%INPUT:
%   Q: joint angles
%   L: link lengths
%OUTPUT:
%   Mij: convertion matrix from j-frame to i-frame
M00_=[cos(pi/4) 0 cos(pi/4) 0
    0 1 0 0
    -cos(pi/4) 0 cos(pi/4) 0
    0 0 0 1];
M0_g=[cos(Q(1)-pi/4) 0 sin(Q(1)-pi/4) 0;
    0 1 0 0;
    -sin(Q(1)-pi/4) 0 cos(Q(1)-pi/4) 0;
    0 0 0 1];
Mg1_=[1 0 0 L(1);
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];
M01_=M00_*M0_g*Mg1_;
M1_g=[cos(-Q(1)+pi/4) 0 sin(-Q(1)+pi/4) 0;
    0 1 0 0;
    -sin(-Q(1)+pi/4) 0 cos(-Q(1)+pi/4) 0;
    0 0 0 1];
Mg1=[1 0 0 L(2)*cos(pi/4);
    0 1 0 0;
    0 0 1 L(2)*cos(pi/4);
    0 0 0 1];
M1_1=M1_g*Mg1;
M11_=[1 0 0 0;
    0 1 0 0;
    0 0 1 L(3)*cos(pi/4);
    0 0 0 1];
M1_2_=[cos(Q(2)) -sin(Q(2)) 0 0;
    sin(Q(2)) cos(Q(2)) 0 0;
    0 0 1 0;
    0 0 0 1];
M12_=M11_*M1_2_;
M2_2=[1 0 0 L(3)*sin(pi/4);
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];
M12=M12_*M2_2;
M23_=[cos(Q(3)-pi/4) 0 sin(Q(3)-pi/4) 0;
        0 1 0 0;
        -sin(Q(3)-pi/4) 0 cos(Q(3)-pi/4) 0;
        0 0 0 1];
M3_3=[1 0 0 L(4);
        0 1 0 0;
        0 0 1 0;
        0 0 0 1];
M23=M23_*M3_3;
end