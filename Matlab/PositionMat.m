function [M01,M12,M23] = PositionMat(Q,L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M01=[1 0 0 L(1)*cos(Q(1))+L(2);
    0 1 0 0;
    0 0 1 -L(1)*sin(Q(1));
    0 0 0 1];
M11_=[cos(pi/4) 0 sin(pi/4) 0;
    0 1 0 0;
    -sin(pi/4) 0 cos(pi/4) 0;
    0 0 0 1];
M1_2_=[cos(Q(2)) -sin(Q(2)) 0 0;
    sin(Q(2)) cos(Q(2)) 0 0;
    0 0 1 0;
    0 0 0 1];
M12_=M11_*M1_2_*inv(M11_);
M2_2=[cos(pi/4) 0 sin(pi/4) L(3)*cos(pi/4);
    0 1 0 0;
    -sin(pi/4) 0 cos(pi/4) L(3)*cos(pi/4);
    0 0 0 1];
M12=M12_*M2_2*inv(M12_);
M23_=[cos(Q(3)) 0 sin(Q(3)) 0;
        0 1 0 0;
        -sin(Q(3)) 0 cos(Q(3)) 0;
        0 0 0 1];
M3_3=[1 0 0 L(4);
        0 1 0 0;
        0 0 1 0;
        0 0 0 1];
M23=M23_*M3_3*inv(M23_);
end

