function  [x,OI,II,time]= PDHG(K,alpha,w_r,tau,sigma,optval,Tol)
    %% PDHG
    x = zeros(size(K,2),1);
    y = zeros(size(K,1),1);
    Kx=K*x;

    OI=0; %%%%%
    II=0; %%%%%
    tic;
while 2>1  %%%%%
        OI=OI+1; %%%%%
        II=II+1; %%%%%
        xpre=x;
        ypre=y;
        Kxpre=Kx;

        x=ProxF(xpre-tau*(K'*ypre),tau,alpha,w_r);
        Kx=K*x;
        y=ProxG(ypre+sigma*(2*Kx-Kxpre));

        
        f = norm(Kx,1)+alpha*(x'*w_r);  %%%%%%
        if -(f-optval)/optval <= Tol  %%%%%%
            time=toc;   %%%%%%%
            break;    %%%%%%
        end   %%%%%%%%
    end
end

