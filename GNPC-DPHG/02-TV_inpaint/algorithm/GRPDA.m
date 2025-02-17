function  [x,SNR,time]= GRPDA(z,I,lamda,psi,A,B,tau,sigma,Tol)
    x=z;
    y=A*x;
    xbar=x;
    i=1;
    
    NI = norm(I(:));
    SNR(i) = 20*log10(NI/(norm(x(:)-I(:)))); 
    time(i)=0;
    i=i+1;
    tic
    while time(i-1) <=1000
        xpre = x;
        ypre = y;
        xbarpre=xbar;
        
        xbar=((psi-1)/psi)*xpre+(1/psi)*xbarpre;

        x=prox_tau_g(xbar - tau * (A' * ypre), tau,lamda,B,z);
        x=min(max(x, 0), 1);

        y=ypre + sigma*(A*x);
        y=min(max(y, -1), 1);

        time(i)=toc;
        SNR(i) = 20*log10(NI/(norm(x(:)-I(:)))); 
        
        Dgap = norm(y(:)-ypre(:))/norm(y(:));
        if Dgap < Tol
            break;  
        end
        i=i+1;
    end
end











