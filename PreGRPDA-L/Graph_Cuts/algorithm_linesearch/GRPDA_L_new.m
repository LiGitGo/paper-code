function  [x,OI,II,time]= GRPDA_L_new(K,alpha,w_r,a,tau,sigma,kesai,optval,Tol)
%%
    x = zeros(size(K,2),1);
    y = zeros(size(K,1),1);
    xbar=x;
    Kty=K'*y;

    OI=0;
    II=0;
    tic;
    while 2>1
        OI=OI+1;
        xpre=x;
        ypre=y;
        xbarpre=xbar;
        taupre=tau;
        sigmapre=sigma;
        Ktypre=Kty;

        xbar=a*xpre+(1-a)*xbarpre;

        x=ProxF(xbar-taupre*(Ktypre),taupre,alpha,w_r);
        Kx=K*x;

        tau=kesai*taupre;
        sigma=kesai*sigmapre;
        while 2>1
            II=II+1;
            %º∆À„y_n
            y=ProxG(ypre+sigma*(Kx));
            
            Kty=K'*y;
            sita=tau/taupre; 
%             if sqrt(tau*sigma)*norm(Kty-Ktypre,2) <=0.99*sqrt(sita/(1-a))*norm((y-ypre),2)
            if sqrt(tau*sigma)*norm(Kty-Ktypre,2) <=0.99*sqrt((2*a*a-6*a+4-sita)*sita/(1-a)/(a*a-3*a+2))*norm((y-ypre),2)
                break
            end
            tau=tau*0.7;
            sigma=sigma*0.7;
        end

        f = norm(Kx,1)+alpha*(x'*w_r);
        if -(f-optval)/optval <= Tol
            time=toc;
            break;
        end
    end
end