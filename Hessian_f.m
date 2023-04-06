function del2F=Hessian_f(xstar,Xref,k) 
    x=xstar;
    n=size(x,1);
    del2F=zeros(n,n);
    h=0.001;
    for i=1:n
        for j=(i+1):n
            e=zeros(n,1);
            e(i)=1;
            e(j)=1;
            A=main_function(x+h*e,Xref,k);
            B=main_function(x-h*e,Xref,k);
            e(i)=-1;
            C=main_function(x+h*e,Xref,k);
            D=main_function(x-h*e,Xref,k);
            del2F(i,j)=(A+B-C-D)/(4*h*h);
            del2F(j,i)=del2F(i,j);
        end
    end        
    for i=1:n
        e=zeros(n,1);
        e(i)=1;
        del2F(i,i)=(main_function(x+h*e,Xref,k)-2*main_function(x,Xref,k)+main_function(x-h*e,Xref,k))/(h*h);
    end
end