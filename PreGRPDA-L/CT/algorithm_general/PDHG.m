function  [x,f,time]= PDHG(A,B,b,tau,sigma)
    %% PDHG
    [rA,cA]=size(A);
    x = zeros(cA,1);
    p = zeros(rA,1);
    q = zeros(size(B,1),1);
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
        Axpre=Ax;
        Bxpre=Bx;

        p=(ppre+sigma*(Axpre-b))/(1+sigma);
        q=Proj(qpre+sigma*Bxpre);
        x=xpre-tau*(A'*(2*p-ppre)+B'*(2*q-qpre));
        Ax=A*x;
        Bx=B*x;

        f(i) = 1/2*sum((Ax-b).^2)+norm(Bx,1);
        time(i)=toc;
        i=i+1;
    end
end
