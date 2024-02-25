%% PDHG
ga=ga_gen;    
step=4/3;  
tau = sqrt(step)/L*ga;
sigma = sqrt(step)/L/ga;
[f_pdhg,t_pdhg] = PDHG(K,b,u,tau,sigma);

%% PDHG-L  
ga=ga_gen;    
step=4/3;  
tau = sqrt(step)/L*ga;
sigma = sqrt(step)/L/ga;
[f_pdhg_L,t_pdhg_L] = PDHG_L(K,b,u,tau,sigma);

%% GRPDA
a=(2-sqrt(2))/2;  
ga=ga_gen;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
tau = sqrt(h_a)/L*ga;
sigma = sqrt(h_a)/L/ga;
[f_grpda,t_grpda]= GRPDA(K,b,u,a,tau,sigma);

%% preGRPDA
a=(2-sqrt(2))/2;  
ga=ga_pre;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
step=h_a;   
gama1=sqrt(1/step)/ga;
gama2=sqrt(1/step)*ga;
M1=gama1*sum(abs(K))';
M2=gama2*sum(abs(K),2);
[f_pregrpda_D,t_pregrpda_D]= pre_GRPDA_D(K,b,u,a,M1,M2);


%% GRPDA-L
a=(2-sqrt(2))/2;  
ga=ga_gen;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
tau = sqrt(h_a)/L*ga;
sigma = sqrt(h_a)/L/ga;
a=1/3;  
[f_grpda_L,t_grpda_L]= GRPDA_L(K,b,u,a,tau,sigma);


%% R-GRPDA-L
a=(2-sqrt(2))/2;  
ga=ga_gen;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
tau = sqrt(h_a)/L*ga;
sigma = sqrt(h_a)/L/ga;
a=0.5; 
kesai=1.1;
[f_grpda_L_new,t_grpda_L_new]= GRPDA_L_new(K,b,u,a,tau,sigma,kesai);

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
[f_pregrpda_L,t_pregrpda_L]= preGRPDA_L(K,b,u,a,M1,M2,kesai);

%% save data
% save('data/L=100.mat')
% save('data/L=300.mat') 
% save('data/L=500.mat') 
% save('data/L=1000.mat') 
% save('data/L=1500.mat') 
% save('data/L=2000.mat') 
