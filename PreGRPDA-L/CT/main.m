%% PDHG
ga=ga_gen;
step=4/3;
tau = sqrt(step)/L*ga; 
sigma = sqrt(step)/L/ga;
[x_pdhg,f_pdhg,t_pdhg]= PDHG(A,B,b,tau,sigma);

%% PDHG-L
ga=ga_gen;
step=4/3;
tau = sqrt(step)/L*ga; 
sigma = sqrt(step)/L/ga;
[x_pdhg_L,f_pdhg_L,t_pdhg_L]= PDHG_L(A,B,b,tau,sigma);

%% GRPDA
a=(2-sqrt(2))/2;  
ga=ga_gen;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
tau = sqrt(h_a)/L*ga;
sigma = sqrt(h_a)/L/ga;
[x_grpda,f_grpda,t_grpda]= GRPDA(A,B,b,a,tau,sigma);

%% ¶Ô½ÇpreGRPDA
a=(2-sqrt(2))/2;  
ga=ga_pre;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
step=h_a;   
gama1=sqrt(1/step)/ga;
gama2=sqrt(1/step)*ga;
M1=gama1*sum(abs(K))';
M2A=gama2*sum(abs(A),2);
M2B=gama2*sum(abs(B),2);
[x_pregrpda,f_pregrpda,t_pregrpda]= preGRPDA(A,B,b,a,M1,M2A,M2B);

%% GRPDA-L-primal
a=(2-sqrt(2))/2;  
ga=ga_gen;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
tau = sqrt(h_a)/L*ga;
sigma = sqrt(h_a)/L/ga;

a=1/3; 
kesai=a*a-3*a+2;
[x_grpda_L_primal,f_grpda_L_primal,t_grpda_L_primal]= GRPDA_L(A,B,b,a,tau,sigma,kesai);

%% GRPDA-L-new
a=(2-sqrt(2))/2;  
ga=ga_gen;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
tau = sqrt(h_a)/L*ga;
sigma = sqrt(h_a)/L/ga;
a=0.5; 
kesai=1.1;
[x_grpda_L,f_grpda_L,t_grpda_L]= GRPDA_L_new(A,B,b,a,tau,sigma,kesai);

%% preGRPDA-L
a=(2-sqrt(2))/2;  
ga=ga_pre;
h_a=(2*a*a-6*a+3)/((1-a)*(a*a-3*a+2));
step=h_a;   
gama1=sqrt(1/step)/ga;
gama2=sqrt(1/step)*ga;
M1=gama1*sum(abs(K))';
M2A=gama2*sum(abs(A),2);
M2B=gama2*sum(abs(B),2);
a=0.5; 
kesai=1.1;
[x_pregrpda_L,f_pregrpda_L,t_pregrpda_L]= preGRPDA_L(A,B,b,a,M1,M2A,M2B,kesai);

%% save data
% save('data/512_50s.mat')