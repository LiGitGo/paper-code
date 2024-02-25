% Proximal operator for the conjudge of g :   
% g(y)=1/2||y-b||^2,   g*(y)=1/2||y||^2+b'y
function  y= prox_sigma_g_conj(v,sigma,b)
    y=1./(1+sigma).*(v-sigma.*b);
end
