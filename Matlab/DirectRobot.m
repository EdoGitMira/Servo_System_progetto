function S = DirectRobot(Q,L)
[M01_ M1_1 M12 M23]=PositionMat(Q,L);
M01=M01_*M1_1;
M02=M01*M12;
M03=M02*M23;
S=M03*[0 0 0 1]';
end

