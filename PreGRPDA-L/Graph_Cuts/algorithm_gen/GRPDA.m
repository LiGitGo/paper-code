function  [x,OI,II,time]= GRPDA(K,alpha,w_r,a,tau,sigma,optval,Tol)
    %% PDHG
    x = zeros(size(K,2),1);
    y = zeros(size(K,1),1);
    xbar=x;

    OI=0; %%%%%
    II=0; %%%%%
    tic;
    while 2>1  %%%%%
        OI=OI+1; %%%%%
        II=II+1; %%%%%
        xpre=x;
        ypre=y;
        xbarpre=xbar;

        xbar=a*xpre+(1-a)*xbarpre;
        x=ProxF(xbar-tau*(K'*ypre),tau,alpha,w_r);
        Kx=K*x;
        y=ProxG(ypre+sigma*Kx);

        f = norm(Kx,1)+alpha*(x'*w_r);  %%%%%%
        if -(f-optval)/optval <= Tol  %%%%%%
            time=toc;   %%%%%%%
            break;    %%%%%%
        end   %%%%%%%%
    end
end