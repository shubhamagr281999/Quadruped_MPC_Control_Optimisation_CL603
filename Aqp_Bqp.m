function [Aqp,Bqp] = Aqp_Bqp(X_ref,n_steps,dt,m,g,I,r,onGround)    
    tempA=eye(13,13);
    Aqp=[];
    Bqp=[];
    for k=1:n_steps
        R=[cos(X_ref(3,k)), sin(X_ref(3,k)), 0;
          -sin(X_ref(3,k)), cos(X_ref(3,k)), 0;
                         0,               0, 1];
        
       I_=(R*I*R')^-1;
                
        A=[zeros(3), zeros(3),        R, zeros(3), zeros(3,1);
           zeros(3), zeros(3), zeros(3),   eye(3), zeros(3,1);
           zeros(3), zeros(3), zeros(3), zeros(3), zeros(3,1);
           zeros(3), zeros(3), zeros(3), zeros(3),          g;
           zeros(1,3),zeros(1,3),zeros(1,3),zeros(1,3),     0];
        B=[];
        UCount=0;
        for leg=1:4
            if(onGround(leg,k)==1)
                B=[B,[ zeros(3);zeros(3);I_*CP(r(:,leg,k));eye(3)./m;zeros(1,3)]];
                UCount=UCount+3;
            end
        end       
       Ai=exp(A.*dt);
       Bi=integral(@(t) exp(A.*t)*B,0,dt,"ArrayValued",true);
       tempA=Ai*tempA;
       Aqp=[Aqp;tempA];
       if(k==1)
           Bqp=Bi;
       else 
           Bqp=[Bqp,zeros((k-1)*13,UCount)];
           Bqp=[Bqp;Ai*Bqp(end-12:end,:)];
           Bqp(end-12:end,(end-UCount+1):end)=Bi;
       end
    end
    function R = CP(r)
        R=[0 -r(3) r(2);r(3) 0 -r(1);-r(2) r(1) 0];
    end
end