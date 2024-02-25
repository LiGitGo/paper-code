function  [x,OI,II,time]= preGRPDA_L(K,alpha,w_r,a,M1,M2,kesai,optval,Tol)
    %% PDHG
    x = zeros(size(K,2),1);
    y = zeros(size(K,1),1);
    xbar=x;
    Kty=K'*y;

    tauM=1./M1;
    sigmaM=1./M2;

    OI=0; 
    II=0; 
    tic;
    while 2>1  
        OI=OI+1; 
        xpre=x;
        ypre=y;
        xbarpre=xbar;
        tauMpre=tauM;
        Ktypre=Kty;

        xbar=a*xpre+(1-a)*xbarpre;
        
        x=ProxF(xbar-tauM.*Ktypre,tauM,alpha,w_r);
        Kx=K*x;
        
        tauM=kesai*tauMpre;
        sigmaM=kesai*sigmaM;
        sitaM=kesai;
        while 2>1
            II=II+1; 
            %º∆À„y_n
            y=ProxG(ypre+sigmaM.*(Kx));
            Kty=K'*y;
            
            if norm(sqrt(tauM).*(Kty-Ktypre),2)<=0.99*sqrt((2*a*a-6*a+4-sitaM)*sitaM/(1-a)/(a*a-3*a+2))*norm(1./sqrt(sigmaM).*(y-ypre),2)
                break;
            end
            tauM=0.7.*tauM;
            sigmaM=0.7*sigmaM;
            sitaM=sitaM*0.7;
        end

        f = norm(Kx,1)+alpha*(x'*w_r);  
        if -(f-optval)/optval <= Tol  
            time=toc;   
            break;    
        end   
    end
end