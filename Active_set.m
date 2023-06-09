function U  = Active_set(x0,Xref,n_steps,mu,fmin,fmax,Q,R,onGround,Bi,Ai)
    %initial guess
    Xref_=reshape(Xref,1,[])';
    xinitial=zeros(24*n_steps,1);
    xPrev=x0;
    for i=0:(n_steps-1)
        for j=0:3
            if(onGround(j+1,i+1)==1)
                xinitial((12*n_steps+i*12+j*3+1):(12*n_steps+i*12+j*3+3))=[0;0;(fmin+fmax)*0.5];
            end
        end
        xinitial((i*12+1):(i*12+12))=Ai(:,:,i+1)*xPrev+Bi(:,:,i+1)*[xinitial((12*(i+n_steps)+1):(12*(i+n_steps)+12));1];
        xPrev=xinitial((i*12+1):(i*12+12));
    end

    onGround=reshape(onGround,1,[])';
    to_keep_constraint=zeros(6*size(onGround,1),1); 
    cons_count=0;
    for i=1:int8(size(onGround,1))
        if(onGround(i)==1)
            to_keep_constraint((cons_count+1):(cons_count+6),1)=((i-1)*6+1):6*i;
            cons_count=cons_count+6;
        end
    end
    
    to_keep_constraints=to_keep_constraint(1:cons_count,1);
    [Du,C,Du_E,C_E] = equality_constraint();
    [I,b_I] = Inequality_cons();
    m_eq=size(Du,1)+size(C,1);
    A = [Du;C;I];
    b = [Du_E;C_E;b_I];
    n = size(A,2);
    m = size(A,1);    
    kk=1; 
    
    x=zeros(24*n_steps,1500);
    p=zeros(24*n_steps,1500);
    xfinal=zeros(24*n_steps,1);
    U=zeros(12,1);
    x(:,kk) = xinitial;
    z = A*x(:,kk)-b;
    W = find(abs(z)<=0.0000001);  
    tW = [1:m];
    while kk<500        
        Anew = [A([W],:)];
        row = size(Anew,1);
        G= Hessian_f(x(:,kk))
        d = grad_f(zeros(size(xinitial)));   % linear part of quadratic function G*xinitial;
        
        kkt = [G Anew';Anew zeros(row,row)]; %check this also
        %kkt = eye(192);
        g = G*x(:,kk) + d;
        Rhs = [g;zeros(row,1)];
        values = inv(kkt)*Rhs;
        p(:,kk) = -values(1:n,:);
        pk = p(:,kk);
        norm_tol = norm(pk);
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
            alphak = min([1,alph])
          
            Wnew = find(alph<=(alphak+10^-3));
            x(:,kk+1) = x(:,kk) + alphak*pk;
            kk = kk+1;
            x(39:3:end,kk)
            W = [W;Wnew'];
        end
    end

    U=zeros(12,1);    
    U(:,1)=x((12*(n_steps)+1):(12*(n_steps)+12),kk-1);

    function [D_,C,D_E,C_E] = equality_constraint()
        C = zeros(12*n_steps,24*n_steps);
        C_E = zeros(12*n_steps,1);
        uConsCount = 1;
        D=zeros(12*n_steps,24*n_steps);
        for j = 1:length(onGround)
            if onGround(j)==0
                D(uConsCount:uConsCount+2,:)=zeros(3,24*n_steps);
                D(uConsCount:uConsCount+2,12*n_steps + 3*j -2: 12*n_steps + 3*j) = eye(3,3);
                uConsCount = uConsCount+3;
            end
        end
        D_=D(1:(uConsCount-1),:);
        D_E=zeros(size(D_,1),1);
        for i=0:(n_steps-1)
            if i==0
                C(1:12,1:12) = eye(12);
                C(1:12,12*(i)+12*n_steps+1:12*(i)+12*n_steps+12) = -Bi(1:12,1:12,i+1);
                C_E(1:12) = Ai(:,:,1)*x0+Bi(:,end,i+1);
            else
                C(12*i+1:12*i+12,12*i+1:12*i+12) = eye(12);
                C(12*i+1:12*i+12,12*(i-1)+1:12*(i-1)+12) = -Ai(1:12,1:12,i+1);
                C(12*i+1:12*i+12,12*i+12*n_steps+1:12*i+12*n_steps+12) = -Bi(1:12,1:12,i+1);
                C_E(12*i+1:12*i+12)=Bi(:,end,i+1);
            end
        end        
    end 

    function [I,b_I] = Inequality_cons()
        I = zeros(24*n_steps,24*n_steps);
        eps=0.01;
        b_I = zeros(24*n_steps,1);
        for i = 0:(n_steps-1)
            for j = 1:4
                if i==0
                    I(6*(j-1)+1,12*i+12*n_steps+1+3*(j-1)+2) = 1;
                    I(6*(j-1)+2,12*i+12*n_steps+1+3*(j-1)+2) = -1;
                    I(6*(j-1)+3,12*i+12*n_steps+1+3*(j-1)+0) = 1 ;
                    I(6*(j-1)+3,12*i+12*n_steps+1+3*(j-1)+2) = mu;
                    I(6*(j-1)+4,12*i+12*n_steps+1+3*(j-1)+0) = -1;
                    I(6*(j-1)+4,12*i+12*n_steps+1+3*(j-1)+2) = mu;
                    I(6*(j-1)+5,12*i+12*n_steps+1+3*(j-1)+1) = 1;
                    I(6*(j-1)+5,12*i+12*n_steps+1+3*(j-1)+2) = mu;
                    I(6*(j-1)+6,12*i+12*n_steps+1+3*(j-1)+1) = -1;
                    I(6*(j-1)+6,12*i+12*n_steps+1+3*(j-1)+2) = mu;
                else 
                    I(24*i+6*(j-1)+1,12*i+12*n_steps+1+3*(j-1)+2) = 1;
                    I(24*i+6*(j-1)+2,12*i+12*n_steps+1+3*(j-1)+2) = -1;
                    I(24*i+6*(j-1)+3,12*i+12*n_steps+1+3*(j-1)+0) = 1 ;
                    I(24*i+6*(j-1)+3,12*i+12*n_steps+1+3*(j-1)+2) = mu;
                    I(24*i+6*(j-1)+4,12*i+12*n_steps+1+3*(j-1)+0) = -1;
                    I(24*i+6*(j-1)+4,12*i+12*n_steps+1+3*(j-1)+2) = mu;
                    I(24*i+6*(j-1)+5,12*i+12*n_steps+1+3*(j-1)+1) = 1;
                    I(24*i+6*(j-1)+5,12*i+12*n_steps+1+3*(j-1)+2) = mu;
                    I(24*i+6*(j-1)+6,12*i+12*n_steps+1+3*(j-1)+1) = -1;
                    I(24*i+6*(j-1)+6,12*i+12*n_steps+1+3*(j-1)+2) = mu;
                end
            end
        end

        for i = 0:n_steps-1
            for j = 1:4
                if i==0
                    b_I(6*(j-1)+1) = fmin;
                    b_I(6*(j-1)+2) = fmax;
                    b_I(6*(j-1)+3) = eps;
                    b_I(6*(j-1)+4) = -eps;
                    b_I(6*(j-1)+5) = eps;
                    b_I(6*(j-1)+6) = -eps;
                else
                    b_I(24*i+6*(j-1)+1) = fmin;
                    b_I(24*i+6*(j-1)+2) = fmax;
                    b_I(24*i+6*(j-1)+3) = eps;
                    b_I(24*i+6*(j-1)+4) = -eps;
                    b_I(24*i+6*(j-1)+5) = eps;
                    b_I(24*i+6*(j-1)+6) = -eps;
                end
            end
        end
        
        I=I(to_keep_constraints,:);
        b_I=b_I(to_keep_constraints,:);
    end

    function del2F=Hessian_f(xstar) 
        x_=xstar;
        n_=size(x_,1);
        del2F=zeros(n_,n_);
        h=0.001;
        for i_=1:int8(n_)
            for j_=(double(i_)+1):double(n_)
                e=zeros(n_,1);
                e(i_)=1;
                e(j_)=1;
                A__=main_function(x_+h*e);
                B__=main_function(x_-h*e);           
                e(i_)=-1;
                C__=main_function(x_+h*e);
                D__=main_function(x_-h*e);
                del2F(i_,j_)=(A__+B__-C__-D__)/(4*h*h);
                del2F(j_,i_)=del2F(i_,j_);
            end
        end        
        for i_=1:n_
            e=zeros(n_,1);
            e(i_)=1;
            del2F(i_,i_)=(main_function(x_+h*e)-2*main_function(x_)+main_function(x_-h*e))/(h*h);
        end
    end

    function model = main_function(X)
        model = 0;
        for i_=0:n_steps-1
            if i_==0
                j_=double(1); 
                p__ = 12*double(n_steps) +1; 
            else
                j_ = double(12*i_+1);
                p__ = double(12*n_steps + 12*(i_) + 1);
            end
            model = model + weighted_norm(X(j_:j_+11) - Xref_(j_:j_+11),Q) + weighted_norm(X(p__:p__+11),R);             
        end
    end

    function N = weighted_norm(a,weight)
        N = a'*weight*a;
    end

    function delF=grad_f(xstar) 
        n_ = length(xstar);     % upto which the definite loop will run for
        h = 0.001;           
        delF=zeros(n_,1);
        for i_ = 1:n_
            e = zeros(n_,1);
            e(i_,1) = 1;     % ei vector for ith component
            xstar_plus_ei_h = xstar+e*h;
            xstar_minus_ei_h = xstar-e*h;
            diff_f_wrt_xi = (main_function(xstar_plus_ei_h)-main_function(xstar_minus_ei_h))/(2*h);    % using central difference
            delF(i_,1) = diff_f_wrt_xi;
        end
    end
end