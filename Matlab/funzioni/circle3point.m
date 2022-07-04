function [X, Xp, Xpp] = circle3point(P1,P2,P3,s,sp,spp)
%INPUT:
%   -P1,P2,P3 tre punti per cui passa la circonferenza
%   -s,sp,spp spostamento, velocità e accelerazione che ha il robot sulla
%   traiettoria, considerando un caso monodimensionale
%OUTPUT:
%   -X,Xp,Xpp variabili dello spazio di lavoro per il punto considerato
%   sulla traiettoria
Xc=centre(P1,P2,P3);
P1(4)=1;
P2(4)=1;
P3(4)=1;
radius=norm(Xc-P1);
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
% if (teta1>teta2 && teta1<teta3) || (teta1<teta2 && teta1<teta3)
%     s=-s;
% end
if max([teta1,teta2,teta3])==teta1
    dteta=teta3-teta1;
else
    dteta=-2*pi+teta3-teta1;
end
tetax=teta1+dteta*s;
Pxc=[radius*cos(tetax);
    radius*sin(tetax);
    0;
    1];
tetaxp=tetax+pi/2*dteta/abs(dteta);
tetaxpp=tetax+pi;
Pxcp=[sp*cos(tetaxp);
    sp*sin(tetaxp);
    0;
    1];
Pxcpp=[spp*cos(tetaxp)+sp^2/radius*cos(tetaxpp);
    spp*sin(tetaxp)+sp^2/radius*sin(tetaxpp);
    0;
    1];
X=M0C*Pxc;
M0C=M0C(:,1:3);
M0C(:,4)=[0 0 0 1]';
Xp=M0C*Pxcp;
Xpp=M0C*Pxcpp;
end
%%
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