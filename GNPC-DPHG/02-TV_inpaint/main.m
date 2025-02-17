%% DPHG
a1=1.1;
b1=0;
a2=2-a1;
b2=0.1;  

step1=4*(1+b1)*(a2+b2)/a1/a2;
step2=(1-b1)*(2-a2-b2)/(1-3/4*a1*a2);
step=min(step1,step2);  

tau= 1;
sigma= step*0.95/tau/L/L;

fprintf('b1 = %.2f -- ¡¾-1,1¡¿\n', b1);
fprintf('a1a2 = %.2f  --  ¡¾0,1.3333¡¿\n', a1*a2);
fprintf('a2+b2 = %.2f -- ¡¾0,2¡¿\n', a2+b2);
fprintf('rs = %.2f\n', 1/tau/sigma);
fprintf('k1 = %.2f\n', 1/step1*L*L);
fprintf('k2 = %.2f\n', 1/step2*L*L);
fprintf('sigma*tau*L*L = %.2f\n', sigma*tau*L*L);

[x_dphg,SNR_dphg,t_dphg] = DPHG(z,I,lamda,A,B,tau,sigma,Tol,a1,a2,b1,b2);

%% PCPDHG
tau = 1;
step = 4;
sigma = step*0.95/(tau*L*L);

[x_pcpdhg,SNR_pcpdhg,t_pcpdhg] = PCPDHG(z,I,lamda,A,B,tau,sigma,Tol);

%% NPCPDHG
a1=1.1;
b1=0;
a2=2-a1;
b2=0.1;    

step1=4*(1+b1)*(a2+b2)/a1/a2;
step2=(1-b1)*(2-a2-b2)/(1-3/4*a1*a2);
step=min(step1,step2);  

tau=1;
sigma= step*0.95/tau/L/L;

fprintf('b1 = %.2f -- ¡¾-1,1¡¿\n', b1);
fprintf('a1a2 = %.2f  --  ¡¾0,1.3333¡¿\n', a1*a2);
fprintf('a2+b2 = %.2f -- ¡¾0,2¡¿\n', a2+b2);
fprintf('rs = %.2f\n', 1/tau/sigma);
fprintf('k1 = %.2f\n', 1/step1*L*L);
fprintf('k2 = %.2f\n', 1/step2*L*L);
fprintf('sigma*tau*L*L = %.2f\n', sigma*tau*L*L);

[x_npcpdhg,SNR_npcpdhg,t_npcpdhg] = NPCPDHG(z,I,lamda,A,B,tau,sigma,Tol,a1,a2,b1,b2);



%% PDHG
tau = 1;   
step =1;
sigma = step*0.95/(tau*L*L);
sita=1;
[x_pdhg,SNR_pdhg,t_pdhg] = PDHG(z,I,lamda,A,B,tau,sigma,Tol,sita);


%% GRPDA

psi = 1.414;
tau = 1;   
step =psi*(2+2*psi-psi*psi)/(1+psi);
sigma = step*0.95/(tau*L*L);
[x_grpda,SNR_grpda,t_grpda]= GRPDA(z,I,lamda,psi,A,B,tau,sigma,Tol);




%% save data
% save('data\result_housergb_256_char_10s_pdhg_1e-4.mat')
% save('data\result_housergb_256_char_10s_pdhg_1e-5.mat')
% save('data\result_fruits_512_row_50s_pdhg_1e-4.mat')
% save('data\result_fruits_512_row_50s_pdhg_1e-5.mat')
