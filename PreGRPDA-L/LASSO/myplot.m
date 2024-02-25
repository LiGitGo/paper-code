%% load data
% load('data/L=100.mat')
% load('data/L=300.mat') 
% load('data/L=500.mat') 
% load('data/L=1000.mat') 
% load('data/L=1500.mat') 
% load('data/L=2000.mat') 
%% plot
figure(1) 
fstar=inf;
fstar = min(fstar,min(f_pdhg));
fstar = min(fstar,min(f_pdhg_L));
fstar = min(fstar,min(f_grpda));
fstar = min(fstar,min(f_pregrpda_D));
fstar = min(fstar,min(f_grpda_L));
fstar = min(fstar,min(f_grpda_L_new));
fstar = min(fstar,min(f_pregrpda_L));

semilogy(t_pdhg,f_pdhg-fstar,'b','DisplayName','PDA','LineWidth', 1)
hold on; legend show;

semilogy(t_pdhg_L,f_pdhg_L-fstar,'g','DisplayName','PDA-L','LineWidth', 1)
hold on; legend show;

semilogy(t_grpda,f_grpda-fstar,'y','DisplayName','GRPDA','LineWidth', 1)
hold on; legend show;

semilogy(t_pregrpda_D,f_pregrpda_D-fstar,'c','DisplayName','PreGRPDA','LineWidth', 1)
hold on; legend show;

semilogy(t_grpda_L,f_grpda_L-fstar,'m','DisplayName','GRPDA-L','LineWidth', 1)
hold on; legend show;

semilogy(t_grpda_L_new,f_grpda_L_new-fstar,'k','DisplayName','R-GRPDA-L','LineWidth', 1)
hold on; legend show;

semilogy(t_pregrpda_L,f_pregrpda_L-fstar,'r','DisplayName','PreGRPDA-L','LineWidth', 1)
hold on; legend show;

xlabel('CPU time, seconds'); ylabel('$\Psi (x^k)-\Psi ^*$','Interpreter','latex');
% xlim([0,10]);

