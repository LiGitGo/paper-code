%% PDHG
ga=ga_gen;
step=4/3;
tau = sqrt(step)/L*ga; 
sigma = sqrt(step)/L/ga;
[x_pdhg,OI_pdhg,II_pdhg,t_pdhg]= PDHG(K,alpha,w_r,tau,sigma,optval,Tol);

%% PDHG-L
ga=ga_gen;
step=4/3;
tau = sqrt(step)/L*ga; 
sigma = sqrt(step)/L/ga;
[x_pdhg_l,OI_pdhg_l,II_pdhg_l,t_pdhg_l]= PDHG_L(K,alpha,w_r,tau,sigma,optval,Tol);


%% GRPDA
a=(2-sqrt(2))/2;  
ga=ga_gen;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
tau = sqrt(h_a)/L*ga;
sigma = sqrt(h_a)/L/ga;
[x_grpda,OI_grpda,II_grpda,t_grpda]= GRPDA(K,alpha,w_r,a,tau,sigma,optval,Tol);


%% ¶Ô½ÇpreGRPDA
a=(2-sqrt(2))/2;  
ga=ga_pre;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
step=h_a;   
gama1=sqrt(1/step)/ga;
gama2=sqrt(1/step)*ga;
M1=gama1*sum(abs(K))';
M2=gama2*sum(abs(K),2);
[x_pregrpda,OI_pregrpda,II_pregrpda,t_pregrpda]= preGRPDA(K,alpha,w_r,a,M1,M2,optval,Tol);

%% GRPDA-L-PRIMAL
a=(2-sqrt(2))/2;  
ga=ga_gen;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
tau = sqrt(h_a)/L*ga;
sigma = sqrt(h_a)/L/ga;

a=1/3;
kesai=a*a-3*a+2;
[x_grpda_L,OI_grpda_L,II_grpda_L,t_grpda_L]= GRPDA_L(K,alpha,w_r,a,tau,sigma,kesai,optval,Tol);


%% GRPDA-L-new
a=(2-sqrt(2))/2;  
ga=ga_gen;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
tau = sqrt(h_a)/L*ga;
sigma = sqrt(h_a)/L/ga;

a=0.5;
kesai=1.1;
[x_grpda_L_new,OI_grpda_L_new,II_grpda_L_new,t_grpda_L_new]= GRPDA_L_new(K,alpha,w_r,a,tau,sigma,kesai,optval,Tol);


%% preGRPDA-L
a=(2-sqrt(2))/2;  
ga=ga_pre;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
step=h_a;   
gama1=sqrt(1/step)/ga;
gama2=sqrt(1/step)*ga;
M1=gama1*sum(abs(K))';
M2=gama2*sum(abs(K),2);

a=0.5;
kesai=1.1;
[x_pregrpda_L,OI_pregrpda_L,II_pregrpda_L,t_pregrpda_L]= preGRPDA_L(K,alpha,w_r,a,M1,M2,kesai,optval,Tol);

% save('data/blueflower_1024_10_1_1e-7_cvx.mat')
% save('data/blueflower_1024_10_1_1e-8_cvx.mat')
% save('data/yellowflower_500_10_1_1e-7_cvx.mat')
% save('data/yellowflower_500_10_1_1e-8_cvx.mat')


