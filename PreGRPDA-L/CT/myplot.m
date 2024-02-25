%% load data
% load('data/512_50s.mat')
%% plot-f
figure(1);  
semilogy(t_pdhg,(f_pdhg-optval)/optval,'b','DisplayName','PDA','LineWidth', 1)
hold on;legend show;

semilogy(t_pdhg_L,(f_pdhg_L-optval)/optval,'g','DisplayName','PDA-L','LineWidth', 1)
hold on;legend show;

semilogy(t_grpda,(f_grpda-optval)/optval,'y','DisplayName','GRPDA','LineWidth', 1)
hold on;legend show;

semilogy(t_pregrpda,(f_pregrpda-optval)/optval,'c','DisplayName','PreGRPDA')
hold on;legend show;

semilogy(t_grpda_L_primal,(f_grpda_L_primal-optval)/optval,'m','DisplayName','GRPDA-L','LineWidth', 1)
hold on;legend show;

semilogy(t_grpda_L,(f_grpda_L-optval)/optval,'k','DisplayName','R-GRPDA-L','LineWidth', 1)
hold on;legend show;

semilogy(t_pregrpda_L,(f_pregrpda_L-optval)/optval,'r','DisplayName','preGRPDA-L')
hold on;legend show;

xlabel('CPU time, seconds');
ylabel('$ {\frac{{\Psi (x^k) - {\Psi ^*}}}{{{\Psi ^*}}}} $','FontSize',13,...
    'Interpreter','latex','FontWeight','bold','FontName','FixedWidth');
%% plot_t
figure(2)
plot_true=reshape(x_true,[M,N]);
imshow(plot_true);

figure(3)
plot_grpda=reshape(x_grpda,[M,N]);
imshow(plot_grpda,[]);

figure(4)
plot_pregrpda_L=reshape(x_pregrpda_L,[M,N]);
imshow(plot_pregrpda_L,[]);

