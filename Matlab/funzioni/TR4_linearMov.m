function [Q,Qp,Qpp,Fq,T,X,Xp,Xpp] = TR4_linearMov(S,law,time,L,M,F,Qi,cons,lambda)
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
    Q=nan;
    Qp=nan;
    Qpp=nan;
    T=nan;
    X=nan;
    Xp=nan;
    Xpp=nan;
    Fq=nan;
    return
end
G=[0 0 -9.81 0 0 -9.81 0 0 -9.81 0 0 -9.81 0 0 -9.81 0 0 -9.81 0 0 0]';
Q=zeros(100,3);
Qp=zeros(100,3);
Qpp=zeros(100,3);
X=zeros(100,4);
Xp=zeros(100,4);
Xpp=zeros(100,4);
Fq=zeros(100,3);
Si=S(:,1);
Sf=S(:,2);
T=linspace(0,time,100);
for t=T
    [x, xp, xpp]=law(t,time,Si(1),Sf(1)-Si(1),lambda);
    [y, yp, ypp]=law(t,time,Si(2),Sf(2)-Si(2),lambda);
    [z, zp, zpp]=law(t,time,Si(3),Sf(3)-Si(3),lambda);
    if ~isWS([x y z]')
        Q=nan;
        Qp=nan;
        Qpp=nan;
        T=nan;
        X=nan;
        Xp=nan;
        Xpp=nan;
        Fq=nan;
        return
    end
    Q(find(T==t),:)=TR4_invNumeric([x y z]',L,1000,cons,Qi);
    if isnan(Q(find(T==t),:))
        Q=nan;
        Qp=nan;
        Qpp=nan;
        T=nan;
        X=nan;
        Xp=nan;
        Xpp=nan;
        Fq=nan;
        return
    end
    q(1)=Q(find(T==t),1);
    q(2)=Q(find(T==t),2);
    q(3)=Q(find(T==t),3);
    Fse=M(q)*G+F;
    Qp(find(T==t),:)=inv(TR4_jacobian(Q(find(T==t),:),L))*[xp;yp;zp];
    qp(1)=Qp(find(T==t),1);
    qp(2)=Qp(find(T==t),2);
    qp(3)=Qp(find(T==t),3);
    Qpp(find(T==t),:)=inv(TR4_jacobian(Q(find(T==t),:),L))*([xpp;ypp;zpp]-TR4_jacobianP(Q(find(T==t),:),Qp(find(T==t),:),L)*qp');
    qpp(1)=Qpp(find(T==t),1);
    qpp(2)=Qpp(find(T==t),2);
    qpp(3)=Qpp(find(T==t),3);
    Je=TR4_Je(q);
    Jep=TR4_Jep(qp,qpp);
    Fq(find(T==t),:)=TR4_Torque(qp',qpp',M(q),Je,Jep,Fse);
    Qi(1)=Q(find(T==t),1);
    Qi(2)=Q(find(T==t),2);
    Qi(3)=Q(find(T==t),3);
    X(find(T==t),:)=[x y z 1]';
    Xp(find(T==t),:)=[xp yp zp 1]';
    Xpp(find(T==t),:)=[xpp ypp zpp 1]';
end
end