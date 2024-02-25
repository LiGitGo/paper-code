%% load data
% load('data/blueflower_1024_10_1_1e-7_cvx.mat')
% load('data/blueflower_1024_10_1_1e-8_cvx.mat')
% load('data/yellowflower_500_10_1_1e-7_cvx.mat')
% load('data/yellowflower_500_10_1_1e-8_cvx.mat')

%% plot result
figure(2)
x=reshape(x_pregrpda_L,[m,n]);
imshow(I+1-x);

%% print
fprintf ('Image = %s,  Tol = %e\n',name,Tol);
fprintf ('PDA:  %d,  %d,  %f\n',OI_pdhg,II_pdhg,t_pdhg);
fprintf ('GRPDA:  %d,  %d,  %f\n',OI_grpda,II_grpda,t_grpda);
fprintf ('PDA-L:  %d,  %d,  %f\n',OI_pdhg_l,II_pdhg_l,t_pdhg_l);
fprintf ('GRPDA-L:  %d,  %d,  %f\n',OI_grpda_L,II_grpda_L,t_grpda_L);
fprintf ('R-GRPDA-L:  %d,  %d,  %f\n',OI_grpda_L_new,II_grpda_L_new,t_grpda_L_new);
fprintf ('PreGRPDA:  %d,  %d,  %f\n',OI_pregrpda,II_pregrpda,t_pregrpda);
fprintf ('PreGRPDA-L:  %d,  %d,  %f\n',OI_pregrpda_L,II_pregrpda_L,t_pregrpda_L);
