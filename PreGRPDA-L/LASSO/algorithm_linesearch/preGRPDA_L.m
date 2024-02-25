function  [f_val,time]= preGRPDA_L(K,b,u,a,M1,M2,kesai)
    x=zeros(size(K,2),1);
    xbar=x;
    y=zeros(size(K,1),1);
   
    Kty=K'*y;

    f_val=[];
    
    tauM=1./M1;
    sigmaM=1./M2;
    
    time=0;
    i=1;
    tic
    while time(1,end)<30
        xpre = x;
        ypre = y;
        xbarpre=xbar;
        tauMpre=tauM;
        Ktypre=Kty;

        %xbar_n
        xbar=a*xpre+(1-a)*xbarpre;
        %x_n
        x=prox_tau_f(xbar - tauM .* Ktypre, tauM, u);
        Kx=K*x;
        
        
        tauM=kesai*tauMpre;
        sigmaM=kesai*sigmaM;
        sitaM=kesai;
        while 2>1
            %y_n
            y=1./(sigmaM+1).*(ypre+sigmaM.*(Kx-b));
            Kty=K'*y;
            
            if norm(sqrt(tauM).*(Kty-Ktypre),2)<=0.99*sqrt((2*a*a-6*a+4-sitaM)*sitaM/(1-a)/(a*a-3*a+2))*norm(1./sqrt(sigmaM).*(y-ypre),2)
                break;
            end
            tauM=0.7.*tauM;
            sigmaM=0.7*sigmaM;
            sitaM=sitaM*0.7;
        end

        f_val(i) = 1/2*(norm(Kx-b ,2)^2) + u*norm(x,1);
        time(i)=toc;
        i=i+1;
    end
end



