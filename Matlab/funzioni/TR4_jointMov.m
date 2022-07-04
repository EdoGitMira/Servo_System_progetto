function [Q_,Qp_,Qpp_,Fq_,T,X_,Xp_,Xpp_] = TR4_jointMov(S,law,time,L,M,F,Qi,cons,lambda)
%INPUT:
%   -S vettore composto da due punti, quello di partenza e arrivo
%   -law legge di moto imposta sulla traiettoria lineare nello spazio di
%   lavoro
%   -time tempo di movimentazione
%   -L lunghezza dei link del robot
%   -M matrice delle masse generalizzate
%   -F forze esterne
%   -Qi terna delle variabili di giunto da utilizzare come partenza per la
%   cinematica inversa
%   -cons vincoli dei giunti
%   -lambda parametro che definisce la percentuale di accelerazione e
%   decelerazione nel caso di legge tre tratti (opzionale nel caso non si
%   usi la legge di moto tre tratti)
%OUTPUT:
%   --Q vettore degli angoli di giunto
%   -Qp vettore delle velocità angolari di giunto
%   -Qp vettore delle accelerazioni angolari di giunto
%   -Fq vettore delle coppie dei giunti
%   -T vettore dei tempi
%   -X vettore delle posizioni assunte nello spazio di lavoro
%   -Xp vettore delle velocità assunte nello spazio di lavoro
%   -Xpp vettore delle accelerazioni assunte nello spazio di lavoro
if nargin<9
    lambda=0;
end
if length(S(1,:))-1 ~= length(time)
    Q_=nan;
    Qp_=nan;
    Qpp_=nan;
    T=nan;
    X_=nan;
    Xp_=nan;
    Xpp_=nan;
    Fq_=nan;
    return
end
G=[0 0 -9.81 0 0 -9.81 0 0 -9.81 0 0 -9.81 0 0 -9.81 0 0 -9.81 0 0 0]';
discretization=100;
Q_in=TR4_invNumeric(S,L,1000,cons);
T=nan(1,discretization);
Q_=nan(discretization,3);
Qp_=nan(discretization,3);
Qpp_=nan(discretization,3);
X_=nan(discretization,4);
Xp_=nan(discretization,4);
Xpp_=nan(discretization,4);
Fq_=nan(discretization,3);
Qf=Q_in(:,2);
dQ=Qf-Qi;
if abs(dQ(1))>abs(dQ(1)+2*pi)
   dQ(1)=dQ(1)+2*pi;
elseif abs(dQ(1))>abs(dQ(1)-2*pi)
    dQ(1)=dQ(1)-2*pi;
end
if abs(dQ(2))>abs(dQ(2)+2*pi)
   dQ(2)=dQ(2)+2*pi;
elseif abs(dQ(2))>abs(dQ(2)-2*pi)
   dQ(2)=dQ(2)-2*pi;
end
if abs(dQ(3))>abs(dQ(3)+2*pi)
    dQ(3)=dQ(3)+2*pi;
elseif abs(dQ(3))>abs(dQ(3)-2*pi)
    dQ(3)=dQ(3)-2*pi;
end
T=linspace(0,time,discretization);
for t=T
    [q1, q1p, q1pp]=law(t,time,Qi(1),dQ(1),lambda);
    [q2, q2p, q2pp]=law(t,time,Qi(2),dQ(2),lambda);
    [q3, q3p, q3pp]=law(t,time,Qi(3),dQ(3),lambda);
    Q_(find(T==t),:)=[q1, q2, q3]';
    Qp_(find(T==t),:)=[q1p, q2p, q3p]';
    Qpp_(find(T==t),:)=[q1pp, q2pp, q3pp]';
    Fse=(M(Q_(find(T==t),:))*G)+F;
    Je=TR4_Je([q1, q2, q3]);
    Jep=TR4_Jep([q1, q2, q3],[q1p, q2p, q3p]);
    Fq_(find(T==t),:)=TR4_Torque([q1p q2p q3p]',[q1pp q2pp q3pp]',M([q1, q2, q3]),Je,Jep,Fse);
    X_(find(T==t),:)=TR4_dirKin([q1 q2 q3]',L);
    Xp_(find(T==t),1:3)=TR4_jacobian([q1, q2, q3],L)*[q1p q2p q3p]';
    Xp_(find(T==t),4)=1;
    Xpp_(find(T==t),1:3)=TR4_jacobianP([q1, q2, q3],[q1p q2p q3p],L)*[q1p q2p q3p]'+TR4_jacobian([q1, q2, q3],L)*[q1pp q2pp q3pp]';
    Xpp_(find(T==t),4)=1;
end
end
%%
function size=real_length(array)
size=1;
max_size=length(array);
while true
    if(isnan(array(size)))
        size=size-1;
        return
    end
    if(size==max_size)
        return
    end
    size=size+1;
end
end