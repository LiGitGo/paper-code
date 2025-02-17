function  [x,f,time]= PCPDHG(K,alpha,w_r,tau,sigma)
    %% PDHG
    x = zeros(size(K,2),1);
    y = zeros(size(K,1),1);
    Kx=K*x;
    
    f = [];
    time=0;
    i=1;
    tic;
    while time(1,end)<300
        xpre=x;
        ypre=y;
        Kxpre=Kx;
       
        xbo=ProxF(xpre-tau*(K'*ypre),tau,alpha,w_r);
        Kxbo = K*xbo;
        ybo=ProxG(ypre+sigma*Kxbo);

        x = xbo+tau/2*(K'*(ypre-ybo));
        Kx = K*x;
        y= ybo +  sigma/2*(2*Kx-Kxpre-Kxbo);
        

        f(i) = norm(Kx,1)+alpha*(x'*w_r);
        time(i)=toc;
        i=i+1;
    end
end

