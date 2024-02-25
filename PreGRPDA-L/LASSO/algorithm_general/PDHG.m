function  [f_val,time]= PDHG(K,b,u,tau,sigma)
    x=zeros(size(K,2),1);
    y=zeros(size(K,1),1);
    Kx=K*x;

    f_val=[];
    time=0;
    i=1;
    tic
    while time(1,end)<30
        xpre = x;
        ypre = y;
        Kxpre=Kx;

        %x_n
        x=prox_tau_f(xpre - tau * (K' * ypre), tau, u);
        Kx=K*x;
        %y_n
        y=prox_sigma_g_conj(ypre+sigma*(2*Kx-Kxpre),sigma,b);

        f_val(i) = 1/2*(norm(Kx-b ,2)^2) + u*norm(x,1);
        time(i)=toc;
        i=i+1;
    end
end












