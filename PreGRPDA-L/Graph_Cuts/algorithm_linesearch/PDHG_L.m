function [x,OI,II,time]= PDHG_L(K,alpha,w_r,tau,sigma,optval,Tol)
    %% pdhg算法求解
    %初始化x0,y0
    x = zeros(size(K,2),1);
    y = zeros(size(K,1),1);
    Kx=K*x;
    Kty=K'*y;
    sita=1;
    
    OI=0; %%%%%
    II=0; %%%%%
    tic;
    while 2>1  %%%%%
        OI=OI+1; %%%%%
        xpre=x;
        ypre=y;
        taupre=tau;
        sigmapre=sigma;
        sitapre=sita;
        Kxpre=Kx;
        Ktypre=Kty;

        x=ProxF(xpre-taupre*(Ktypre),taupre,alpha,w_r);
        Kx=K*x;

        tau=taupre*sqrt(1+sitapre);
        sigma=sigmapre*sqrt(1+sitapre);
        while 2>1
            II=II+1; %%%%%%%%
            sita=tau/taupre;
            Kxbar=(1+sita)*Kx-sita*Kxpre;
            %计算y_n
            y=ProxG(ypre+sigma*(Kxbar));
            
            Kty=K'*y;
            if sqrt(tau*sigma)*norm(Kty-Ktypre,2) <=0.99*norm((y-ypre),2)
                break
            end
            tau=tau*0.7;
            sigma=sigma*0.7;
        end

        f = norm(Kx,1)+alpha*(x'*w_r);  %%%%%%
        if -(f-optval)/optval <= Tol  %%%%%%
            time=toc;   %%%%%%%
            break;    %%%%%%
        end   %%%%%%%%
    end
end










