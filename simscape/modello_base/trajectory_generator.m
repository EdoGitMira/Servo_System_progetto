function [output] = trajectory_generator(t)
q1=0;
q2=1;
q3=2;
q1p=3;
q2p=4;
q3p=5;
q1pp=6;
q2pp=7;
q3pp=8;
Fq1=9;
Fq2=10;
Fq3=11;
output =[q1,q2,q3,q1p,q2p,q3p,q1pp,q2pp,q3pp,Fq1,Fq2,Fq3];
end

