function [I,b_I] = Inequality_cons(k,mu,fmin,fmax,pos)
    I = zeros(24*k,24*k);
    b_I = zeros(24*k,1);
        for i = 0:k-1
            for j = 1:4
                if i==0
                    I(6*(j-1)+1,12*i+12*k+1+3*(j-1)+2) = 1;
                    I(6*(j-1)+2,12*i+12*k+1+3*(j-1)+2) = -1;
                    I(6*(j-1)+3,12*i+12*k+1+3*(j-1)+0) = 1 ;
                    I(6*(j-1)+3,12*i+12*k+1+3*(j-1)+2) = mu;
                    I(6*(j-1)+4,12*i+12*k+1+3*(j-1)+0) = -1;
                    I(6*(j-1)+4,12*i+12*k+1+3*(j-1)+2) = mu;
                    I(6*(j-1)+5,12*i+12*k+1+3*(j-1)+1) = 1;
                    I(6*(j-1)+5,12*i+12*k+1+3*(j-1)+2) = mu;
                    I(6*(j-1)+6,12*i+12*k+1+3*(j-1)+1) = -1;
                    I(6*(j-1)+6,12*i+12*k+1+3*(j-1)+2) = mu;
                else 
                    I(24*i+6*(j-1)+1,12*i+12*k+1+3*(j-1)+2) = 1;
                    I(24*i+6*(j-1)+2,12*i+12*k+1+3*(j-1)+2) = -1;
                    I(24*i+6*(j-1)+3,12*i+12*k+1+3*(j-1)+0) = 1 ;
                    I(24*i+6*(j-1)+3,12*i+12*k+1+3*(j-1)+2) = mu;
                    I(24*i+6*(j-1)+4,12*i+12*k+1+3*(j-1)+0) = -1;
                    I(24*i+6*(j-1)+4,12*i+12*k+1+3*(j-1)+2) = mu;
                    I(24*i+6*(j-1)+5,12*i+12*k+1+3*(j-1)+1) = 1;
                    I(24*i+6*(j-1)+5,12*i+12*k+1+3*(j-1)+2) = mu;
                    I(24*i+6*(j-1)+6,12*i+12*k+1+3*(j-1)+1) = -1;
                    I(24*i+6*(j-1)+6,12*i+12*k+1+3*(j-1)+2) = mu;
                end
            end
        end

    for i = 0:k-1
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
    to_keep_constraints=[];
    for i=1:size(pos,1)
        if(pos(i)==1)
            to_keep_constraints=[to_keep_constraints,((i-1)*6+1):6*i];
        end
    end
    I=I(to_keep_constraints,:);
    b_I=b_I(to_keep_constraints,:);
end