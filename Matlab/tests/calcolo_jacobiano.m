%% Calcolo Jacobiano
clear all
%definizione variabili simboliche
Q=sym("q",[1 3]).'; 
L=sym("l",[1 4]);
%calcolo posizione rispetto alle variabili di giunto e alle lunghezze dei
%link
S=TR4_dirKin(Q,L)
%calcolo jacobiana
J=jacobian(S,Q)
J=matlabFunction(J)
load TR4_data.mat
S=TR4_dirKin(Q,link)
J=jacobian(S,Q)
syms Sq1 Sq2 Sq3 Cq1 Cq2 Cq3
J=J(1:3,:)
J=subs(J,sin(Q(1)),Sq1);
J=subs(J,sin(Q(2)),Sq2);
J=subs(J,sin(Q(3)),Sq3);
J=subs(J,cos(Q(1)),Cq1);
J=subs(J,cos(Q(2)),Cq2);
J=subs(J,cos(Q(3)),Cq3)
latex(J)