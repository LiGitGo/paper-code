function  [x,SNR,time]= PDHG(z,I,lamda,A,B,tau,sigma,Tol,sita)
    x=z;
    y=A*x;
    i=1;
    
    NI = norm(I(:));
    SNR(i) = 20*log10(NI/(norm(x(:)-I(:)))); 
    time(i)=0;
    i=i+1;
    tic
    while time(i-1) <=1000
        xpre = x;
        ypre = y;
        
        x=prox_tau_g(xpre - tau * (A' * ypre), tau,lamda,B,z);
        x=min(max(x, 0), 1);

        xbar=x+sita*(x-xpre);

        y=ypre + sigma*(A*xbar);
        y=min(max(y,-1), 1);

        time(i)=toc;
        SNR(i) = 20*log10(NI/(norm(x(:)-I(:)))); 
        
        Dgap = norm(y(:)-ypre(:))/norm(y(:));
        if Dgap < Tol
            break;  
        end
        i=i+1;
    end
end











