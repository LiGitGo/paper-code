function  [x,f,time]= GRPDA(K,alpha,w_r,psi,tau,sigma)
%% PDHG
x = zeros(size(K,2),1);
y = zeros(size(K,1),1);
xbar=x;

f = [];
time=0;
i=1;
tic;
while time(1,end)<300
    xpre=x;
    ypre=y;
    xbarpre=xbar;
    
    xbar=((psi-1)/psi)*xpre+(1/psi)*xbarpre;
    x=ProxF(xbar-tau*(K'*ypre),tau,alpha,w_r);
    Kx=K*x;
    y=ProxG(ypre+sigma*Kx);
    
    f(i) = norm(Kx,1)+alpha*(x'*w_r);
    time(i)=toc;
    i=i+1;
end