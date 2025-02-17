function  [x,f,time]=NPCPDHG(K,alpha,w_r,tau,sigma,a1,a2,b1,b2)
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

        x = (1+b1)*xbo-b1*xpre+tau*a1/2*(K'*(ypre-ybo));
        Kx = K*x;
        y= ypre -(a2+b2)*(ypre-ybo) +  sigma*a2/2*(2*Kx-Kxpre-Kxbo);
        

        f(i) = norm(Kx,1)+alpha*(x'*w_r);
        time(i)=toc;
        i=i+1;
    end
end

