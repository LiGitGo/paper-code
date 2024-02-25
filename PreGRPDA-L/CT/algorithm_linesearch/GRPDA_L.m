function  [x,f,time]= GRPDA_L(A,B,b,a,tau,sigma,kesai)
    %% PDHG
    [rA,cA]=size(A);
    x = zeros(cA,1);
    p = zeros(rA,1);
    q = zeros(size(B,1),1);
    pbar=p;
    qbar=q;
    Ax=A*x;
    Bx=B*x;

    f = []; 
    time=0;
    i=1;
    tic;
    while time(1,end)<50
        xpre=x;
        ppre=p;
        qpre=q;
        pbarpre=pbar;
        qbarpre=qbar;
        taupre=tau;
        sigmapre=sigma;
        Axpre=Ax;
        Bxpre=Bx;

        pbar=a*ppre+(1-a)*pbarpre;
        qbar=a*qpre+(1-a)*qbarpre;
        p=(pbar+sigmapre*(Axpre-b))/(1+sigmapre);
        q=Proj(qbar+sigmapre*(Bxpre));
        Atp=A'*p;
        Btq=B'*q;
        AAtp=A*Atp;
        ABtq=A*Btq;
        BAtp=B*Atp;
        BBtq=B*Btq;

        tau=kesai*taupre;
        sigma=kesai*sigmapre;
        while 2>1
            x=xpre-tau*(Atp+Btq);
            Ax=Axpre-tau*(AAtp+ABtq);
            Bx=Bxpre-tau*(BAtp+BBtq);

            sita=tau/taupre;
            if sqrt(tau*sigma)*norm([Ax;Bx]-[Axpre;Bxpre],2) <=0.99*sqrt(sita/(1-a))*norm((x-xpre),2)
                break
            end
            tau=tau*0.7;
            sigma=sigma*0.7;
        end

        f(i) = 1/2*sum((Ax-b).^2)+norm(Bx,1);
        time(i)=toc;
        i=i+1;
    end
end
