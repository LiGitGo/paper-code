function  [x,f,time]= preGRPDA_L(A,B,b,a,M1,M2A,M2B,kesai)
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
        tauMpre=tauM;
        sigmaMApre=sigmaMA;
        sigmaMBpre=sigmaMB;
        Axpre=Ax;
        Bxpre=Bx;

        pbar=a*ppre+(1-a)*pbarpre;
        qbar=a*qpre+(1-a)*qbarpre;

        p=(pbar+sigmaMApre.*(Axpre-b))./(1+sigmaMApre);
        q=Proj(qbar+sigmaMBpre.*(Bxpre));
        Atp=A'*p;
        Btq=B'*q;
        MpreAtp=tauMpre.*(Atp+Btq);
        AMAtp=A*MpreAtp;
        BMAtp=B*MpreAtp;

        tauM=kesai*tauMpre;
        sigmaMA=kesai*sigmaMA;
        sigmaMB=kesai*sigmaMB;
        sitaM=kesai;
        while 2>1
            x=xpre-tauM.*(Atp+Btq);
            Ax=Axpre-sitaM*AMAtp;
            Bx=Bxpre-sitaM*BMAtp;
            
            if norm(sqrt([sigmaMA;sigmaMB]).*([Ax;Bx]-[Axpre;Bxpre]),2)<=0.99*sqrt((2*a*a-6*a+4-sitaM)*sitaM/(1-a)/(a*a-3*a+2))*norm(1./sqrt(tauM).*(x-xpre),2)
                break
            end

            tauM=tauM*0.7;
            sigmaMA=sigmaMA*0.7;
            sigmaMB=sigmaMB*0.7;
            sitaM=sitaM*0.7;
        end

        f(i) = 1/2*sum((Ax-b).^2)+norm(Bx,1);
        time(i)=toc;
        i=i+1;
    end
end