function [Q,Qp,Qpp,Fq,T,X,Xp,Xpp]=TR4_move(S,T,law,L,M,F,mode,cons,Qi,lambda)
%INPUT:
%   -S punti della movimentazione
%   -T tempo di movimentazione (non utilizzato nel caso di movimentazione
%   con tempo di attuazione minimo)
%   -law funzione della legge oraria da seguire
%   -L lunghezza dei link
%   -M matrice generalizzata delle masse
%   -F forze esterne
%   -mode tipologia di movimentazione:
%       *linear movimentazione lineare tra due punti
%       *joint movimentazione nello spazio dei giunti
%       *circular movimentazione circolare con un punto di partenza uno di
%       arrivo e uno per definire la circonferenza desiderata
%       *minimum_time_TS movimentazione nello spazio dei giunti con tempo
%       minimo di movimentazione considerando la legge di moto tre tratti e
%       le limitazioni di accelerazione e velocità dei giunti
%       *minimum_time_C movimentazione nello spazio dei giunti con tempo
%       minimo di movimentazione considerando la legge di moto cicloidale e
%       le limitazioni di accelerazione e velocità dei giunti
%   -cons vincoli dei giunti
%   -Qi terna di variabili dei giunti (opzionale, se non fornita viene
%   cercata all'interno dello spazio dei giunti raggiungibile)
%   -lambda parametro che definisce la percentuale di accelerazione e
%   decelerazione nel caso di legge tre tratti (opzionale nel caso non si
%   usi la legge di moto tre tratti)
%OUTPUT:
%   -Q vettore degli angoli di giunto
%   -Qp vettore delle velocità angolari di giunto
%   -Qp vettore delle accelerazioni angolari di giunto
%   -Fq vettore delle coppie dei giunti
%   -T vettore dei tempi
%   -X vettore delle posizioni assunte nello spazio di lavoro
%   -Xp vettore delle velocità assunte nello spazio di lavoro
%   -Xpp vettore delle accelerazioni assunte nello spazio di lavoro
if nargin<10
    lambda=0;
end
if  ~isWS(S)
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
switch mode
    case 'linear'
        if nargin<9
            Qi=TR4_invNumeric(S(:,1),L,1000,cons);
        end
        [Q,Qp,Qpp,Fq,T,X,Xp,Xpp]=TR4_linearMov(S,law,T,L,M,F,Qi,cons,lambda);
    case 'joint'
        if nargin<9
            Qi=TR4_invNumeric(S(:,1),L,1000,cons);
        end
        [Q,Qp,Qpp,Fq,T,X,Xp,Xpp]=TR4_jointMov(S,law,T,L,M,F,Qi,cons,lambda);
    case 'circular'
        if nargin<9
            Qi=TR4_invNumeric(S(:,1),L,1000,cons);
        end
        [Q,Qp,Qpp,Fq,T,X,Xp,Xpp]=TR4_circularMov(S,law,T,L,M,F,Qi,cons,lambda);
    case 'minimum_time_TS'
        if length(S(1,:))-1 ~= length(T)
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
        if nargin<9
            Qi=TR4_invNumeric(S(:,1),L,1000,cons);
        end
        Q_in=Qi;
        Q_in=[Q_in TR4_invNumeric(S(2,:),L,1000,cons)];
        dQ=Q_in(:,1)-Q_in(:,2);
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
        [T_min1,l1] = minimumTimeTS(abs(dQ(1)),cons(2,1),cons(3,1));
        [T_min2,l2] = minimumTimeTS(abs(dQ(2)),cons(2,2),cons(3,2));
        [T_min3,l3] = minimumTimeTS(abs(dQ(3)),cons(2,3),cons(3,3));
        T_min=max([T_min1,T_min2,T_min3]);
        if T_min1>=T_min
            l=l1;
        elseif T_min2>=T_min
            l=l2;
        else
            l=l3;
        end
        [Q,Qp,Qpp,Fq,T,X,Xp,Xpp]=TR4_jointMov(S,@tretratti,T_min,L,M,F,Qi,cons,l);
    case 'minimum_time_C'
        if length(S(1,:))-1 ~= length(T)
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
        if nargin<9
            Qi=TR4_invNumeric(S(:,1),L,1000,cons);
            if isnan(Qi)
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
        end
        Q_in=Qi;
        Q_in=[Q_in TR4_invNumeric(S(:,2),L,1000,cons)];
        dQ=Q_in(:,1)-Q_in(:,2);
        T_min1 = minimumTimeC(abs(dQ(1)),cons(2,1),cons(3,1));
        T_min2 = minimumTimeC(abs(dQ(2)),cons(2,2),cons(3,2));
        T_min3 = minimumTimeC(abs(dQ(3)),cons(2,3),cons(3,3));
        T_min=max([T_min1,T_min2,T_min3]);
        [Q,Qp,Qpp,Fq,T,X,Xp,Xpp]=TR4_jointMov(S,@Cicloidale,T_min,L,M,F,Qi,cons);
    otherwise
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
end
%%
function [T_min,l] = minimumTimeTS(dq,vmax,amax)
t1=sqrt(dq/amax);
v1=t1*amax;
if v1>vmax
    v1=vmax;
end
t1=v1/amax;
s1=amax*0.5*t1^2;
if s1*2>dq
    s1=dq/2;
    t1=sqrt(2*s1/amax);
    v1=t1*amax;
end
t2=(dq-2*s1)/v1;
T_min=t1*2+t2;
l=t1/T_min;
end
function T_min = minimumTimeC(dq,vmax,amax)
Tv=2*dq/vmax;
Ta=sqrt(2*pi*dq/amax);
T_min=max([Tv,Ta]);
end