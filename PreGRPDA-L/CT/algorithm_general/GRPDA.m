function  [x,f,time]= GRPDA(A,B,b,a,tau,sigma)
    %% GRPDA
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
        pbarpre=pbar;
        qbarpre=qbar;
        ppre=p;
        qpre=q;
        Axpre=Ax;
        Bxpre=Bx;

        pbar=a*ppre+(1-a)*pbarpre;
        qbar=a*qpre+(1-a)*qbarpre;
        p=(pbar+sigma*(Axpre-b))/(1+sigma);
        q=Proj(qbar+sigma*Bxpre);
        x=xpre-tau*(A'*p+B'*q);
        Ax=A*x;
        Bx=B*x;


        f(i) = 1/2*sum((Ax-b).^2)+norm(Bx,1);
        time(i)=toc;
        i=i+1;
    end
end