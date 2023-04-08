function [xfinal,x] = Active_set(n_steps,Xref,onGround,mu,fmin,fmax,x0,Q,R,Ai,Bi)

    %initial guess
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

    to_keep_constraints=[]; 
    for i=1:size(onGround,1)
        if(onGround(i)==1)
            to_keep_constraints=[to_keep_constraints,((i-1)*6+1):6*i];
        end
    end

    [D,C,D_E,C_E] = equality_constraint();
    [I,b_I] = Inequality_cons();
    A = [D;C;I];
    b = [D_E;C_E;b_I];
    n = size(A,2);
    m = size(A,1);
    G= Hessian_f(xinitial);
    d = grad_f(zeros(size(xinitial)));   % linear part of quadratic function G*xinitial;
    kk=1;
    x(:,kk) = xinitial;

    z = A*x(:,kk)-b;
    W = find(abs(z)<=0.0000001);
    tW = [1:m];

    while kk<1500
        Anew = [A([W],:)];
        row = size(Anew,1);
        kkt = [G Anew';Anew zeros(row,row)];
        %kkt = eye(192);
        g = G*x(:,kk) + d;
        Rhs = [g;zeros(row,1)];
        values = inv(kkt)*Rhs;
        p(:,kk) = -values(1:n,:);
        pk = p(:,kk);
        norm_tol = norm(pk);
        if norm(pk)<0.0001 
            lambda = values(n+1:length(values),:);
            if all(lambda>=0)
                xfinal = x(:,kk);
                break
            else 
                j=find(lambda==min(lambda));
                W(j) = [];
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
            alph(alph<=0) = 10;
            alphak = min([1,alph]);

            Wnew = find(alph==alphak);
            x(:,kk+1) = x(:,kk) + alphak*pk;
            kk = kk+1;
            W = [W;Wnew'];
        end
    end

    function [D,C,D_E,C_E] = equality_constraint()
        C = zeros(12*n_steps,24*n_steps);
        C_E = zeros(12*n_steps,1);
        p = 1;
        for j = 1:length(onGround)
            if onGround(j)==0
                D(p:p+2,:)=zeros(3,24*n_steps);
                D(p:p+2,12*n_steps + 3*j -2: 12*n_steps + 3*j) = eye(3,3);
                p = p+3;
            end
        end
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
        D_E=zeros(size(D,1),1);
    end 

    function [I,b_I] = Inequality_cons()
        I = zeros(24*n_steps,24*n_steps);
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
                else
                    b_I(24*i+6*(j-1)+1) = fmin;
                    b_I(24*i+6*(j-1)+2) = fmax;
                end
            end
        end
        
        I=I(to_keep_constraints,:);
        b_I=b_I(to_keep_constraints,:);
    end

    function del2F=Hessian_f(xstar) 
        x=xstar;
        n_=size(x,1);
        del2F=zeros(n_,n_);
        h=0.001;
        for i=1:n
            for j=(i+1):n_
                e=zeros(n_,1);
                e(i)=1;
                e(j)=1;
                A_=main_function(x+h*e);
                B_=main_function(x-h*e);
                e(i)=-1;
                C_=main_function(x+h*e);
                D_=main_function(x-h*e);
                del2F(i,j)=(A_+B_-C_-D_)/(4*h*h);
                del2F(j,i)=del2F(i,j);
            end
        end        
        for i=1:n
            e=zeros(n,1);
            e(i)=1;
            del2F(i,i)=(main_function(x+h*e)-2*main_function(x)+main_function(x-h*e))/(h*h);
        end
    end

    function model = main_function(X)
        model = 0;
        for i=0:n_steps-1
            if i==0
                j=1; 
                p = 12*n_steps +1; 
            else
                j = 12*i+1;
                p = 12*n_steps + 12*(i) + 1;
            end
            model = model + weighted_norm(X(j:j+11) - Xref(j:j+11),Q) + weighted_norm(X(p:p+11),R);             
        end
    end

    function N = weighted_norm(a,weight)
        N = a'*weight*a;
    end

    function delF=grad_f(xstar) 
        n_ = length(xstar);     % upto which the definite loop will run for
        h = 0.001;           
        delF=zeros(n_,1);
        for i = 1:n_
            e = zeros(n,1);
            e(i,1) = 1;     % ei vector for ith component
            xstar_plus_ei_h = xstar+e*h;
            xstar_minus_ei_h = xstar-e*h;
            diff_f_wrt_xi = (main_function(xstar_plus_ei_h)-main_function(xstar_minus_ei_h))/(2*h);    % using central difference
            delF(i,1) = diff_f_wrt_xi;
        end
    end

end
