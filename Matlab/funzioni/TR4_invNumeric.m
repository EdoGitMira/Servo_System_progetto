function Q = TR4_invNumeric(S,L,Nmax,cons,Qi)
%INPUT:
%   -S posizzione nello spazio di lavoro di cui bisogna calcolare le
%   variabili di giunto
%   -L lunghezze dei link del robot
%   -Nmax numero massimo di iterazioni
%   -cons limiti in spostamento
%   -Qi terna di variabili di giunto iniziale. 
%   Opzionale: nel caso non venga fornita la funzione cerca la terna
%   desiderata provando 8000 terne di partenze all'interno dello spazio 
%   dei giunti definito dai limiti
%OUTPUT:
%   -Q terna delle variabili di giunto che produce la posizione richiesta,
%   nel caso non sia stata trovata restituisce NaN
Q=nan(3,length(S(1,:)));
index=1;
if nargin<5
    for s=S
        for q1 = linspace(-cons(1,1)/2,cons(1,1)/2,20)
            for q2 = linspace(-cons(1,2)/2,cons(1,2)/2,20)
                for q3 = linspace(-cons(1,3)/2,cons(1,3)/2,20)
                    Qi=[q1 q2 q3]';
                    for i=1:Nmax/10
                        Si=TR4_dirKin(Qi,L);
                        Si=Si(1:end-1);
                        J=TR4_jacobian(Qi,L);
                        if abs(det(J)) <= eps
                            break
                        end
                        Qi=Qi+inv(J)*(s-Si);
                        err=s-Si;
                        Qi(1)=mod(Qi(1)+pi,2*pi)-pi;
                        Qi(2)=mod(Qi(2)+pi,2*pi)-pi;
                        Qi(3)=mod(Qi(3)+pi,2*pi)-pi;
                        if norm(err)<= 10^-10 && Qi(1)>-cons(1,1)/2 && Qi(1)<cons(1,1)/2 && Qi(2)>-cons(1,2)/2 && Qi(2)<cons(1,2)/2 && Qi(3)>-cons(1,3)/2 && Qi(3)<cons(1,3)/2
                            Q(1,index)=mod(Qi(1)+pi,2*pi)-pi;
                            Q(2,index)=mod(Qi(2)+pi,2*pi)-pi;
                            Q(3,index)=mod(Qi(3)+pi,2*pi)-pi;
                            break
                        elseif norm(err)<= 10^-10
                            break
                        end
                    end
                    if norm(err)<= 10^-10 && Qi(1)>-cons(1,1)/2 && Qi(1)<cons(1,1)/2 && Qi(2)>-cons(1,2)/2 && Qi(2)<cons(1,2)/2 && Qi(3)>-cons(1,3)/2 && Qi(3)<cons(1,3)/2
                       break 
                    end
                end
                if norm(err)<= 10^-10 && Qi(1)>-cons(1,1)/2 && Qi(1)<cons(1,1)/2 && Qi(2)>-cons(1,2)/2 && Qi(2)<cons(1,2)/2 && Qi(3)>-cons(1,3)/2 && Qi(3)<cons(1,3)/2
                   break 
                end
            end
            if norm(err)<= 10^-10 && Qi(1)>-cons(1,1)/2 && Qi(1)<cons(1,1)/2 && Qi(2)>-cons(1,2)/2 && Qi(2)<cons(1,2)/2 && Qi(3)>-cons(1,3)/2 && Qi(3)<cons(1,3)/2
               break 
            end
        end
        if norm(err)> 10^-10 || ~(Qi(1)>-cons(1,1)/2 && Qi(1)<cons(1,1)/2 && Qi(2)>-cons(1,2)/2 && Qi(2)<cons(1,2)/2 && Qi(3)>-cons(1,3)/2 && Qi(3)<cons(1,3)/2)
            return
        end
        index=index+1;
    end
else
    for s=S
        for i=1:Nmax
            Si=TR4_dirKin(Qi,L);
            Si=Si(1:end-1);
            J=TR4_jacobian(Qi,L);
            if abs(det(J)) <= eps
                break
            end
            Qi=Qi+inv(J)*(s-Si);
            err=s-Si;
            if norm(err)<= 10^-10 && Qi(1)>-cons(1,1)/2 && Qi(1)<cons(1,1)/2 && Qi(2)>-cons(1,2)/2 && Qi(2)<cons(1,2)/2 && Qi(3)>-cons(1,3)/2 && Qi(3)<cons(1,3)/2
                Q(1,index)=mod(Qi(1)+pi,2*pi)-pi;
                Q(2,index)=mod(Qi(2)+pi,2*pi)-pi;
                Q(3,index)=mod(Qi(3)+pi,2*pi)-pi;
                break 
            end
        end
        if norm(err)> 10^-10
            Q=nan;
            return
        end
        index=index+1;
    end
end
end