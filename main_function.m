function model = main_function(X)
    model = 0;
    for i=0:k-1
        if i==0
            j=1; 
            p = 12*k +1; 
        else
            j = 12*i+1;
            p = 12*k + 12*(i) + 1;
        end
        model = model + weighted_norm(X(j:j+11) - Xref(j:j+11)) + weighted_norm(X(p:p+11));             
    end
end

