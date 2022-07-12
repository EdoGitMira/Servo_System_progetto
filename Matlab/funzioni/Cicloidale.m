function [x,xp,xpp]=Cicloidale(t,T,S0,dS,l1)
%legge di moto Cicloidale (acc.costante)
%
%
% t tempo per cui calcolare la legge
% T tempo di azionamento
% S0 posizione iniziale
% dS ampiezza movimento
% 
% si assume Vini=Vfin=0
%
x=S0+dS/T*(t-T/2/pi*sin(2*pi/T*t));
xp=dS/T*(1-cos(2*pi/T*t));
xpp=dS/T*2*pi/T*sin(2*pi/T*t);
end

