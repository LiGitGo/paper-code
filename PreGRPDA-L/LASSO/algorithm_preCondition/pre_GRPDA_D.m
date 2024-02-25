function  [f_val,time]= pre_GRPDA_D(K,b,u,a,M1,M2);
    x=zeros(size(K,2),1);
    xbar=x;
    y=zeros(size(K,1),1);

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

        %xbar_n
        xbar=a*xpre+(1-a)*xbarpre;
        %x_n
        x=prox_tau_f(xbar - tauM .* (K' * ypre), tauM, u);
        Kx=K*x;
        %y_n
        y=prox_sigma_g_conj(ypre+sigmaM.*Kx,sigmaM,b);

        f_val(i) = 1/2*(norm(Kx-b ,2)^2) + u*norm(x,1);
        time(i)=toc;
        i=i+1;
    end
end











