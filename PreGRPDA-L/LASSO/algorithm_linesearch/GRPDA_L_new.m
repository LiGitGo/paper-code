function  [f_val,time]= GRPDA_L_new(K,b,u,a,tau,sigma,kesai)
    x=zeros(size(K,2),1);
    xbar=x;
    y=zeros(size(K,1),1);
    Kty=K'*y;
    ktb=K'*b;

    f_val=[];
    
    time=0;
    i=1;
    tic
    while time(1,end)<30
        xpre = x;
        ypre = y;
        xbarpre=xbar;
        taupre=tau;
        sigmapre=sigma;
        Ktypre=Kty;
        
        %xbar_n
        xbar=a*xpre+(1-a)*xbarpre;
        
        %x_n
        x=prox_tau_f(xbar - taupre*(Ktypre), taupre, u);
        Kx=K*x;
        ktkx=K'*Kx;
        
        tau=kesai*taupre;
        sigma=kesai*sigmapre;
        while 2>1
            %y_n
            y=(ypre+sigma*(Kx-b))/(1+sigma);
            
            sita=tau/taupre;
            Kty=(Ktypre+sigma*(ktkx-ktb))/(1+sigma);
            if sqrt(tau*sigma)*norm(Kty-Ktypre,2) <=0.99*sqrt((2*a*a-6*a+4-sita)*sita/(1-a)/(a*a-3*a+2))*norm((y-ypre),2)
                break
            end
            tau=tau*0.7;
            sigma=sigma*0.7;
        end
        
        f_val(i) = 1/2*(norm(Kx-b ,2)^2) + u*norm(x,1);
        time(i)=toc;
        i=i+1;
    end
end















