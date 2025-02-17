%% GNPC-DPHG
%算法参数
a1=1.1;
b1=0;
a2=2-a1;
b2=0.35;    

step1=4*(1+b1)*(a2+b2)/a1/a2;
step2=(1-b1)*(2-a2-b2)/(1-3/4*a1*a2);
step=min(step1,step2)*0.95;  

tau = 1;
sigma = step/tau/L/L;


fprintf('b1 = %.2f\n', b1);
fprintf('a1a2 = %.2f\n', a1*a2);
fprintf('a2+b2 = %.2f\n', a2+b2);
fprintf('rs = %.2f\n', 1/tau/sigma);
fprintf('k1 = %.2f\n', 1/step1*L*L);
fprintf('k2 = %.2f\n', 1/step2*L*L);
fprintf('sigma*tau*L*L = %.2f\n', sigma*tau*L*L);
[x_dphg,f_dphg,t_dphg]= DPHG(K,alpha,w_r,tau,sigma,a1,a2,b1,b2);


%% CP-PDA
tau = 1;
step = 1*0.95;
sigma = step/tau/L/L;


[x_pdhg,f_pdhg,t_pdhg]= PDHG(K,alpha,w_r,tau,sigma);


%% GRPDA
psi = 1.414;
tau = 1;   
step =psi*(2+2*psi-psi*psi)/(1+psi);
sigma = step*0.95/(tau*L*L);
[x_grpda,f_grpda,t_grpda]= GRPDA(K,alpha,w_r,psi,tau,sigma);


%% PC-PDHG
tau =1;
step = 4*0.95;
sigma = step/tau/L/L;
[x_pcpdhg,f_pcpdhg,t_pcpdhg]= PCPDHG(K,alpha,w_r,tau,sigma);

%% NPC-PDHG
%算法参数
a1=1.1;
b1=0;
a2=2-a1;
b2=0.35;    

step1=4*(1+b1)*(a2+b2)/a1/a2;
step2=(1-b1)*(2-a2-b2)/(1-3/4*a1*a2);
step=min(step1,step2)*0.95; 

tau = 1;
sigma = step/tau/L/L;

fprintf('b1 = %.2f\n', b1);
fprintf('a1a2 = %.2f\n', a1*a2);
fprintf('a2+b2 = %.2f\n', a2+b2);
fprintf('rs = %.2f\n', 1/tau/sigma);
fprintf('k1 = %.2f\n', 1/step1*L*L);
fprintf('k2 = %.2f\n', 1/step2*L*L);
fprintf('sigma*tau*L*L = %.2f\n', sigma*tau*L*L);
[x_npcpdhg,f_npcpdhg,t_npcpdhg]= NPCPDHG(K,alpha,w_r,tau,sigma,a1,a2,b1,b2);




%% save data
% save('data/yellow_flower_beta10_alpha1_tau1_500s_pdhg_1.mat');
% save('data/blue_flower_300_beta10_alpha1_tau1_500s_pdhg_1.mat');