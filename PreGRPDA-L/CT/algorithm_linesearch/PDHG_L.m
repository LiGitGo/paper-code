function   [x,f,time]= PDHG_L(A,B,b,tau,sigma)
    [rA,cA]=size(A);
    x = zeros(cA,1);
    p = zeros(rA,1);
    q = zeros(size(B,1),1);
    Ax=A*x;
    Bx=B*x;
    Atp=A'*p;
    Btq=B'*q;
    AAtp=A*Atp;
    ABtq=A*Btq;
    BAtp=B*Atp;
    BBtq=B*Btq;
    sita=1;
    
    f = []; 
    time=0;
    i=1;
    tic;
    while time(1,end)<50
        xpre=x;
        ppre=p;
        qpre=q;
        taupre=tau;
        sigmapre=sigma;
        sitapre=sita;
        Axpre=Ax;
        Bxpre=Bx;
        Atppre=Atp;
        Btqpre=Btq;
        AAtppre=AAtp;
        ABtqpre=ABtq;
        BAtppre=BAtp;
        BBtqpre=BBtq;
        

        p=(ppre+sigmapre*(Axpre-b))/(1+sigmapre);
        q=Proj(qpre+sigmapre*(Bxpre));
        Atp=A'*p;
        Btq=B'*q;
        AAtp=A*Atp;
        ABtq=A*Btq;
        BAtp=B*Atp;
        BBtq=B*Btq;
        
        tau=taupre*sqrt(1+sitapre);
        sigma=sigmapre*sqrt(1+sitapre);
        while 2>1
            sita=tau/taupre;
            
            x=xpre-tau*((1+sita)*Atp-sita*Atppre+(1+sita)*Btq-sita*Btqpre);
            
            Ax=Axpre-tau*((1+sita)*AAtp-sita*AAtppre+(1+sita)*ABtq-sita*ABtqpre);
            Bx=Bxpre-tau*Bxpre-tau*((1+sita)*BAtp-sita*BAtppre+(1+sita)*BBtq-sita*BBtqpre);

            if sqrt(tau*sigma)*norm([Ax;Bx]-[Axpre;Bxpre],2) <=0.99*norm((x-xpre),2)
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












