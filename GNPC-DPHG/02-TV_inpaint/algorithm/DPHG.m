function  [x,SNR,time]= DPHG(z,I,lamda,A,B,tau,sigma,Tol,a1,a2,b1,b2)
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
        
        ybo = ypre + tau*(A*xpre);
        ybo=min(max(ybo,-1), 1);

        xbo = prox_tau_g(xpre - sigma * (A' * ybo), sigma,lamda,B,z);
        xbo=min(max(xbo, 0), 1);

        y=(1+b1)*ybo - b1*ypre - (tau*a1)/2*(A*(xpre-xbo));
        x = xpre - (a2+b2)*(xpre-xbo)-(sigma*a2)/2*(A'*(2*y-ypre-ybo)); 
        
        time(i)=toc;
        SNR(i) = 20*log10(NI/(norm(x(:)-I(:)))); 
        
        Dgap = norm(y(:)-ypre(:))/norm(y(:));
        if Dgap < Tol
            break;  
        end
        i=i+1;
    end
end











