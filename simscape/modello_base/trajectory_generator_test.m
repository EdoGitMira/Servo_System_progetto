function [output] = trajectory_generator_test(t)
    if t(1)<T(1) 
       output =[0,0,0,0,0,0,0,0,0,1,0,0];
    elseif t(1)>T(1) && t<T(end-2)
        
        q1=Ql(1,t(2),1);
        q2=Ql(1,t(2),2);
        q3=Ql(1,t(2),3);

        q1p=Ql1p(t(2));
        q2p=Ql2p(t(2));
        q3p=Ql3p(t(2));

        q1pp=Ql1p(t(t(2)));
        q2pp=Ql2p(t(t(2)));
        q3pp=Ql3p(t(t(2)));
        
        if t(1)>T(t(2)) && t(1)>T(t(2)+1)
            t(2) = t(2)+1;
        end
        output =[q1,q2,q3,q1p,q2p,q3p,q1pp,q2pp,q3pp,t(2),0,0]; 
    else
        q1=Ql(1,t(end),1);
        q2=Ql(1,t(end),2);
        q3=Ql(1,t(end),3);

        q1p=Ql1p(t(end));
        q2p=Ql2p(t(end));
        q3p=Ql3p(t(end));

        q1pp=Ql1p(t(end));
        q2pp=Ql2p(t(end));
        q3pp=Ql3p(t(end));
        output =[q1,q2,q3,q1p,q2p,q3p,q1pp,q2pp,q3pp,0,0,0];    
    end
end