function  [x,SNR,time]= NPCPDHG(z,I,lamda,A,B,tau,sigma,Tol,a1,a2,b1,b2)
    
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
        
        x =(1+b1)* xbo - b1*xpre+tau*a1/2*(A'*(ypre-ybo));
        y= ypre -(a2+b2)* (ypre-ybo) + sigma*a2/2*(A*(2*x-xbo-xpre));

        time(i)=toc;
        SNR(i) = 20*log10(NI/(norm(x(:)-I(:)))); 
        
        Dgap = norm(y(:)-ypre(:))/norm(y(:));
        if Dgap < Tol
            break;  
        end
        i=i+1;
    end
end


