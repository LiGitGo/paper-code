function  [x,f,time]= DPHG(K,alpha,w_r,tau,sigma,a1,a2,b1,b2)
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
        Kxpre = Kx;

        ybo=ProxG(ypre+tau*Kxpre);
        xbo=ProxF(xpre-sigma*(K'*ybo),sigma,alpha,w_r);

        y=(1+b1)*ybo - b1*ypre - (tau*a1)/2*(K*(xpre-xbo));

        x = xpre - (a2+b2)*(xpre-xbo)-(sigma*a2)/2*(K'*(2*y-ypre-ybo));
        Kx = K*x;

        f(i) = norm(Kx,1)+alpha*(x'*w_r);
        time(i)=toc;
        i=i+1;
    end
end

