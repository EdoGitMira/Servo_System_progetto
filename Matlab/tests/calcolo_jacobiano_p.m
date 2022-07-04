clear all
%definizione delle variabili di giunto dipendenti dal tempo
syms t q1(t) q2(t) q3(t)
Q(1)=q1;
Q(2)=q2;
Q(3)=q3;
Q=Q.';
%definizione delle variabili che descrivono le lunghezze dei giunti
L=sym("l",[1 4]);
%calcolo della matrice jacobiana attraverso la funzione TR4_jacobian()
J=TR4_jacobian(Q,L)
%differenziazione rispetto al tempo della matrice appena calcolata
Jp=diff(J,t)
syms Q1 Q2 Q3 Q1p Q2p Q3p
Jp=subs(Jp,diff(q2),Q2p);
Jp=subs(Jp,diff(q1),Q1p);
Jp=subs(Jp,diff(q3),Q3p);
Jp=subs(Jp,q1,Q1);
Jp=subs(Jp,q2,Q2);
Jp=subs(Jp,q3,Q3)
Jp=matlabFunction(Jp)
load TR4_data.mat
JL=TR4_jacobian(Q,link)
JpL=diff(JL,t)
syms Sq1 Sq2 Sq3 Cq1 Cq2 Cq3
JpL=subs(JpL,sin(Q(1)),Sq1);
JpL=subs(JpL,sin(Q(2)),Sq2);
JpL=subs(JpL,sin(Q(3)),Sq3);
JpL=subs(JpL,cos(Q(1)),Cq1);
JpL=subs(JpL,cos(Q(2)),Cq2);
JpL=subs(JpL,cos(Q(3)),Cq3)
T=sym('teta',[1 13])
t1=(Sq2*Sq3*diff(Q(3), t))/50;
JpL=subs(JpL,t1,T(1));
t2=(Sq2*Sq3*diff(Q(2), t))/50;
JpL=subs(JpL,t2,T(2));
t3=(Cq3*Sq2*diff(Q(3), t))/50;
JpL=subs(JpL,t3,T(3));
t4=(Cq2*Sq3*diff(Q(3), t))/50;
JpL=subs(JpL,t4,T(4));
t5=(Cq3*Sq2*diff(Q(2), t))/50;
JpL=subs(JpL,t5,T(5));
t6=(Cq2*Sq3*diff(Q(2), t))/50;
JpL=subs(JpL,t6,T(6));
t7=(Cq2*Cq3*diff(Q(3), t))/50;
JpL=subs(JpL,t7,T(7));
t8=(Cq2*Cq3*diff(Q(2), t))/50;
JpL=subs(JpL,t8,T(8));
t9=(399*Cq2*diff(Q(2), t))/25000;
JpL=subs(JpL,t9,T(9));
t10=(Sq3*diff(Q(3), t))/50;
JpL=subs(JpL,t10,T(10));
t11=(Cq3*diff(Q(3), t))/50
JpL=subs(JpL,t11,T(11));
t12=sym(pi)/4 + Q(3);
JpL=subs(JpL,t12,T(12))
t13=(Cq2*cos(T(12))*diff(q3(t), t))/25 - (Sq2*((2*sin(T(12)))/25 + sym('6505604140384725/144115188075855872'))*diff(q2(t), t))/2;
JpL=subs(JpL,t13,T(13))
latex(JpL)