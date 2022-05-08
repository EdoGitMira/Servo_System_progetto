function  plotRobot(Q,L,fig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[M01_ M1_1 M12 M23]=PositionMat(Q,L);
M01=M01_*M1_1;
M02=M01*M12;
M03=M02*M23;
P0=[0 0 0 1];
P1_=M01_*P0';
P1=M01*P0';
P2=M02*P0';
P3=M03*P0';
figure(fig);
plot3([P0(1) P1_(1) P1(1) P2(1) P3(1)],[P0(2) P1_(2) P1(2) P2(2) P3(2)],[P0(3) P1_(3) P1(3) P2(3) P3(3)])
xlabel("X")
ylabel("Y")
zlabel("Z")
axis equal
end

