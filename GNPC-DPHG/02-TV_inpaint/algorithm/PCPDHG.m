function  [x,SNR,time]= PCPDHG(z,I,lamda,A,B,tau,sigma,Tol)
    x=z;
    y=A*x;
    i=1;
    
    NI = norm(I(:));
    SNR(i) = 20*log10(NI/(norm(x(:)-I(:)))); 
    time(i)=0;
    i=i+1;
    tic
    while time(i-1) <= 1000
        xpre = x;
        ypre = y;
        
        xbo=prox_tau_g(xpre - tau * (A' * ypre), tau,lamda,B,z);
        xbo=min(max(xbo, 0), 1);

        ybo=ypre + sigma*(A*xbo);
        ybo=min(max(ybo,-1), 1);
        
        x = xbo+tau/2*(A'*(ypre-ybo));
        y= ybo + sigma/2*(A*(2*x-xbo-xpre));

        time(i)=toc;
        SNR(i) = 20*log10(NI/(norm(x(:)-I(:)))); 
        
        Dgap = norm(y(:)-ypre(:))/norm(y(:));
        if Dgap < Tol
            break;  
        end
        i=i+1;
    end
end


