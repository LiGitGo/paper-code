function  [f_val,time]= PDHG_L(K,b,u,tau,sigma)
    x=zeros(size(K,2),1);
    y=zeros(size(K,1),1);
    kty=K'*y;
    kx=K*x;
    ktb=K'*b;
    ktkx=K'*kx;
    
    f_val=[];
    sita=1;
    
    time=0;
    i=1;
    tic
    while time(1,end)<30
        xpre = x;
        ypre = y;
        taupre=tau;
        sigmapre=sigma;
        sitapre=sita;
        ktypre=kty;
        kxpre=kx;
        ktkxpre=ktkx;

        %x_n
        x=prox_tau_f(xpre - taupre * ktypre, taupre, u);
        kx=K*x;
        ktkx=K'*kx;
        
        tau=taupre*sqrt(1+sitapre);
        sigma=sigmapre*sqrt(1+sitapre);
        while 2>1
            sita=tau/taupre;
            %xbar_n
            xbar = x+sita*(x - xpre);
            Kxbar=(1+sita)*kx-sita*kxpre;
            ktkxbar=(1+sita)*ktkx-sita*ktkxpre;
            %y_n
            y=(ypre+sigma*(Kxbar-b))/(1+sigma);
            kty=(ktypre+sigma*(ktkxbar-ktb))/(1+sigma);
            
            if sqrt(tau*sigma)*norm(kty-ktypre,2) <=0.99*norm((y-ypre),2)
                break
            end
            tau=tau*0.7;
            sigma=sigma*0.7;
        end
        

        f_val(i) = 1/2*(norm(kx-b ,2)^2) + u*norm(x,1);
        time(i)=toc;
        i=i+1;
    end
end












