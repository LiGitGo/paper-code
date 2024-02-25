function  [x,f,time]= preGRPDA(A,B,b,a,M1,M2A,M2B)
    [rA,cA]=size(A);
    x = zeros(cA,1);
    p = zeros(rA,1);
    q = zeros(size(B,1),1);
    pbar=p;
    qbar=q;
    Ax=A*x;
    Bx=B*x;

    tauM=1./M1;
    sigmaMA=1./M2A;
    sigmaMB=1./M2B;

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
        p=(pbar+sigmaMA.*(Axpre-b))./(1+sigmaMA);
        q=Proj(qbar+sigmaMB.*Bxpre);
        x=xpre-tauM.*(A'*p+B'*q);
        Ax=A*x;
        Bx=B*x;


        f(i) = 1/2*sum((Ax-b).^2)+norm(Bx,1);
        time(i)=toc;
        i=i+1;
    end
end
