function delF=grad_f(xstar,Xref,k) 
    n = length(xstar);     % upto which the definite loop will run for
    h = 0.001;           

    for i = 1:n
        e = zeros(n,1);
        e(i,1) = 1;     % ei vector for ith component
        xstar_plus_ei_h = xstar+e*h;
        xstar_minus_ei_h = xstar-e*h;
        diff_f_wrt_xi = (main_function(xstar_plus_ei_h,Xref,k)-main_function(xstar_minus_ei_h,Xref,k))/(2*h);    % using central difference
        delF(i,1) = diff_f_wrt_xi;
    end
end


