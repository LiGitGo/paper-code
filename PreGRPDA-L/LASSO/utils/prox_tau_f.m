% Proximal operator for L1_norm : f=u||x||_1
function  x= prox_tau_f(v,tau,u)
    x=max(0, abs(v)-tau*u).*sign(v);
end