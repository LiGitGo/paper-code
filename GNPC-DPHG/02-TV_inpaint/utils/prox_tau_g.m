function  y = prox_tau_g(v,tau,lamda,B,z)
    y= (tau*lamda*B.*z+v)./(tau*lamda.*B.*B+ones(size(z,1),size(z,2)));
end