function [Q,Qp,Qpp,Fq,T,X,Xp,Xpp] = TR4_circularMov(S,law,time,L,M,F,Qi,cons,lambda)
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
if length(S(1,:)) ~= 3 || length(time) ~= 1
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

   P1=S(:,1);
   P2=S(:,2);
   P3=S(:,3);
   T=linspace(0,time,100);
   disp=circleDisp(S);
   for t=T
        [displ, displp, displpp]=law(t,time,0,disp,lambda);
        fff=displ/disp;
        [Xc, Xpc, Xppc] = circle3point(P1,P2,P3,fff,displp,displpp);
        Q(find(T==t),:)=TR4_invNumeric(Xc(1:3),L,1000,cons,Qi);
        q(1)=Q(find(T==t),1);
        q(2)=Q(find(T==t),2);
        q(3)=Q(find(T==t),3);
        Fse=M(q)*G+F;
        Qp(find(T==t),:)=inv(TR4_jacobian(Q(find(T==t),:),L))*Xpc(1:3);
        qp(1)=Qp(find(T==t),1);
        qp(2)=Qp(find(T==t),2);
        qp(3)=Qp(find(T==t),3);
        Qpp(find(T==t),:)=inv(TR4_jacobian(Q(find(T==t),:),L))*(Xppc(1:3)-TR4_jacobianP(Q(find(T==t),:),Qp(find(T==t),:),L)*qp');
        qpp(1)=Qpp(find(T==t),1);
        qpp(2)=Qpp(find(T==t),2);
        qpp(3)=Qpp(find(T==t),3);
        Je=TR4_Je(q);
        Jep=TR4_Jep(qp,qpp);
        Fq(find(T==t),:)=TR4_Torque(qp',qpp',M(q),Je,Jep,Fse);
        Qi(1)=Q(find(T==t),1);
        Qi(2)=Q(find(T==t),2);
        Qi(3)=Q(find(T==t),3);
        X(find(T==t),:)=Xc;
        Xp(find(T==t),:)=Xpc;
        Xpp(find(T==t),:)=Xppc;
        S=TR4_dirKin(Q(find(T==t),:)',L);
   end

end
%%
function displ=circleDisp(S)
P1=S(:,1);
P2=S(:,2);
P3=S(:,3);
P1(4)=1;
P2(4)=1;
P3(4)=1;
Xc=centre(P1,P2,P3);
radius=norm(Xc-P1);
displ=radius*pi*2;
M0C=transpositionMAT(P1,P2,P3,Xc);
MC0=inv(M0C);
P1c=MC0*P1;
P2c=MC0*P2;
P3c=MC0*P3;
teta1=atan2(P1c(2),P1c(1));
if teta1<0
    teta1=teta1+2*pi;
end
teta2=atan2(P2c(2),P2c(1));
if teta2<0
    teta2=teta2+2*pi;
end
teta3=atan2(P3c(2),P3c(1));
if teta3<0
    teta3=teta3+2*pi;
end
if max([teta1,teta2,teta3])==teta1
    dteta=teta3-teta1;
else
    dteta=-2*pi+teta3-teta1;
end
displ=displ*abs(dteta)/(2*pi);
end

function rx = Xrotation(P1,P2,P3,Pc)

x1=P1(1);
x2=P2(1);
x3=P3(1);
y1=P1(2);
y2=P2(2);
y3=P3(2);
z1=P1(3);
z2=P2(3);
z3=P3(3);
xc=Pc(1);
yc=Pc(2);
zc=Pc(3);
Z=-(y1.*z2-y2.*z1-y1.*z3+y3.*z1+y2.*z3-y3.*z2-x1.*y2.*z3+x1.*y3.*z2+x2.*y1.*z3-x2.*y3.*z1-x3.*y1.*z2+x3.*y2.*z1+x1.*y2.*zc-x1.*yc.*z2-x2.*y1.*zc+x2.*yc.*z1+xc.*y1.*z2-xc.*y2.*z1-x1.*y3.*zc+x1.*yc.*z3+x3.*y1.*zc-x3.*yc.*z1-xc.*y1.*z3+xc.*y3.*z1+x2.*y3.*zc-x2.*yc.*z3-x3.*y2.*zc+x3.*yc.*z2+xc.*y2.*z3-xc.*y3.*z2)./(x1.*y2-x2.*y1-x1.*y3+x3.*y1+x2.*y3-x3.*y2);
rx=atan2(-Z,1);

end

function ry = Yrotation(P1,P2,P3,Pc,rx)

x1=P1(1);
x2=P2(1);
x3=P3(1);
y1=P1(2);
y2=P2(2);
y3=P3(2);
z1=P1(3);
z2=P2(3);
z3=P3(3);
xc=Pc(1);
yc=Pc(2);
zc=Pc(3);
Z=(x1.*z2-x2.*z1-x1.*z3+x3.*z1+x2.*z3-x3.*z2+x1.*y2.*z3-x1.*y3.*z2-x2.*y1.*z3+x2.*y3.*z1+x3.*y1.*z2-x3.*y2.*z1-x1.*y2.*zc+x1.*yc.*z2+x2.*y1.*zc-x2.*yc.*z1-xc.*y1.*z2+xc.*y2.*z1+x1.*y3.*zc-x1.*yc.*z3-x3.*y1.*zc+x3.*yc.*z1+xc.*y1.*z3-xc.*y3.*z1-x2.*y3.*zc+x2.*yc.*z3+x3.*y2.*zc-x3.*yc.*z2-xc.*y2.*z3+xc.*y3.*z2)./(x1.*y2.*cos(rx)-x2.*y1.*cos(rx)-x1.*y3.*cos(rx)+x3.*y1.*cos(rx)+x2.*y3.*cos(rx)-x3.*y2.*cos(rx)+y1.*z2.*sin(rx)-y2.*z1.*sin(rx)-y1.*z3.*sin(rx)+y3.*z1.*sin(rx)+y2.*z3.*sin(rx)-y3.*z2.*sin(rx));
ry=atan2(Z,1);

end

function M0C = transpositionMAT(P1,P2,P3,Pc)
rx=Xrotation(P1,P2,P3,Pc);
ry = Yrotation(P1,P2,P3,Pc,rx);
xc=Pc(1);
yc=Pc(2);
zc=Pc(3);
M0C=reshape([cos(rx),0.0,-sin(rx),0.0,sin(rx).*sin(ry),cos(ry),cos(rx).*sin(ry),0.0,cos(ry).*sin(rx),-sin(ry),cos(rx).*cos(ry),0.0,xc,yc,zc,1.0],[4,4]);
end