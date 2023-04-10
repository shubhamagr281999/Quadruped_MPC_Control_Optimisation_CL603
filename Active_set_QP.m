function [U,Xf]  = Active_set_QP(x0,Xref,n_steps,mu,fmin,fmax,Q_,R_,onGround,Bqp,Aqp)
    %initial guess
    Xref_=reshape(Xref,1,[])';
    onGround=reshape(onGround,1,[])';
    R=zeros(size(Bqp,2),size(Bqp,2));
    legOnGroundCount=size(Bqp,2)/3;
    for i=1:legOnGroundCount
        R(((i-1)*3+1):3*i,((i-1)*3+1):3*i)=R_;
    end
    Q=zeros(13*n_steps,13*n_steps);
    for i=1:n_steps
        Q(((i-1)*13+1):13*i,((i-1)*13+1):13*i)=Q_;
    end
    xinitial=zeros(size(Bqp,2),1);
    xinitial(3:3:end)=(fmin+fmax)*0.5;
    m_eq=0;
    [I,b_I] = Inequality_cons();
    A = I;
    b = b_I;
    n = size(A,2);
    m = size(A,1);    
    kk=1; 
    
    G=(Bqp'*Q*Bqp+R).*2;
    d=Bqp'*Q*(Aqp*x0-Xref_).*2;
    x(:,kk) = xinitial;
    z = A*x(:,kk)-b;
    W = find(abs(z)<=0.0000001);  
    tW = [1:m];    
   
    while kk<500        
        Anew = [A([W],:)];
        row = size(Anew,1);     
        kkt = [G Anew';Anew zeros(row,row)];
        g = G*x(:,kk) + d;
        Rhs = [g;zeros(row,1)];
        values = inv(kkt)*Rhs;
        p(:,kk) = -values(1:n,:);
        pk = p(:,kk);
        if norm(pk)<0.001 
            lambda = values(n+1:length(values),:);
            lambda_=lambda((m_eq+1):size(lambda,1));
            if all(lambda_>=0) 
                break
            else 
                ConsToRemove=find(lambda==min(lambda));
                for i=ConsToRemove
                    if(W(ConsToRemove)>m_eq)
                        W(ConsToRemove) = [];
                    end
                end
                x(:,kk+1) = x(:,kk);
                kk = kk+1;
            end
        else
            nw =  setdiff(tW,W);
            alph = [];
            for i=nw
                if A(i,:)*pk < 0
                    alph(i) = (b(i) - A(i,:)*x(:,kk))/(A(i,:)*pk);               
                end
            end
            alph(alph<=0.00001) = 10;
            alphak = min([1,alph]);
          
            Wnew = find(alph<=(alphak+10^-3));
            x(:,kk+1) = x(:,kk) + alphak*pk;
            kk = kk+1;
            W = [W;Wnew'];
        end
    end
    xfinal=x(:,kk);
    U=zeros(12,1);
    count=1;
    for i=1:4
        if(onGround(i)==1)
            U(((i-1)*3+1):i*3)=xfinal(count:count+2);
            count=count+3;
        end
    end
    Xf=x;    

    function [I,b_I] = Inequality_cons()
        I = zeros(2*size(Bqp,2),size(Bqp,2));
        b_I = zeros(size(I,1),1);
        legOnGroundCount=size(Bqp,2)/3;
        eps=0.01;        
        for i = 1:legOnGroundCount    
            I(6*(i-1)+1,(i-1)*3+3) = 1;
            I(6*(i-1)+2,(i-1)*3+3) = -1;
            I(6*(i-1)+3,(i-1)*3+1) = 1 ;
            I(6*(i-1)+3,(i-1)*3+3) = mu;
            I(6*(i-1)+4,(i-1)*3+1) = -1;
            I(6*(i-1)+4,(i-1)*3+3) = mu;
            I(6*(i-1)+5,(i-1)*3+2) = 1;
            I(6*(i-1)+5,(i-1)*3+3) = mu;
            I(6*(i-1)+6,(i-1)*3+2) = -1;
            I(6*(i-1)+6,(i-1)*3+3) = mu;
            b_I(6*(i-1)+1) = fmin;
            b_I(6*(i-1)+2) = fmax;
            b_I(6*(i-1)+3) = eps;
            b_I(6*(i-1)+4) = -eps;
            b_I(6*(i-1)+5) = eps;
            b_I(6*(i-1)+6) = -eps;
                
        end
    end
end