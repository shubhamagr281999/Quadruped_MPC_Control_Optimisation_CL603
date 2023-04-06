function [D,C,D_E,C_E] = equality_constraint(k,pos,x0)
    D = zeros(12*k,24*k); % on-groun matrice
    D_E = zeros(12*k,1);
    C = zeros(12*k,24*k);
    C_E = zeros(12*k,1);
    p = 1;
        for j = 1:length(pos)
            if pos(j)==0
                D(p:p+2,12*k + 3*j -2: 12*k + 3*j) = eye(3,3);
                p = p+3;
            end
        end
    for i=0:k-1
        if i==0
            C(1:12,1:12) = eye(12);
            C(1:12,12*(i)+12*k+1:12*(i)+12*k+12) = -Bi(1:12,1:12,i+1);
        else
            C(12*i+1:12*i+12,12*i+1:12*i+12) = eye(12);
            C(12*i+1:12*i+12,12*(i-1)+1:12*(i-1)+12) = -Ai(1:12,1:12,i+1);
            C(12*i+1:12*i+12,12*i+12*k+1:12*i+12*k+12) = -Bi(1:12,1:12,i+1);
        end
    end
    C_E(1:12) = A(:,:,1)*x0;
end 
