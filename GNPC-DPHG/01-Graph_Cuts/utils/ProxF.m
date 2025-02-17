function  [x]= ProxF(v,tau,alpha,w_r)
    x=min(max(v-alpha*tau.*w_r,0),1);
end