%% data
% load('data\result_housergb_256_char_10s_pdhg_1e-4.mat')
% load('data\result_housergb_256_char_10s_pdhg_1e-5.mat')
% load('data\result_fruits_512_row_50s_pdhg_1e-4.mat')
% load('data\result_fruits_512_row_50s_pdhg_1e-5.mat')
%% plot-SNR-t
figure(3)

plot(t_dphg,SNR_dphg,'r','DisplayName','GNPC-DPHG')
hold on
legend show;

plot(t_pcpdhg,SNR_pcpdhg,'b','DisplayName','PC-PDHG')
hold on
legend show;

plot(t_npcpdhg,SNR_npcpdhg,'y','DisplayName','NPC-PDHG')
hold on
legend show;


plot(t_pdhg,SNR_pdhg,'g','DisplayName','CP-PDA')
hold on
legend show;

plot(t_grpda,SNR_grpda,'k','DisplayName','GRPDA')
hold on
legend show;

xlabel('Time, s'); ylabel('SNR');
if strcmp(ImageName,'housergb.png')  
    xlim([0,2]);
elseif strcmp(ImageName,'fruits.png')
    xlim([0,40]);
end
legend( 'Location','southeast');



%% plot-SNR-it
figure(4)

plot(SNR_dphg,'r','DisplayName','GNPC-DPHG')
hold on
legend show;

plot(SNR_pcpdhg,'b','DisplayName','PC-PDHG')
hold on
legend show;

plot(SNR_npcpdhg,'y','DisplayName','NPC-PDHG')
hold on
legend show;

plot(SNR_pdhg,'g','DisplayName','CP-PDA')
hold on
legend show;

plot(SNR_grpda,'k','DisplayName','GRPDA')
hold on
legend show;

xlabel('Iteration, No.'); ylabel('SNR');
if strcmp(ImageName,'housergb.png')  
    xlim([0,200]);
elseif strcmp(ImageName,'fruits.png')
    xlim([0,800]);
end
legend( 'Location','southeast');

%% plot-result
figure(5)
imshow(reshape(x_dphg, [n, n, 3]));

%% print
fprintf ('Image=%s,  Tol=%e\n',ImageName,Tol);
fprintf ('PDHG:  %d/%f/%f\n',size(t_pdhg,2)-1,t_pdhg(1,end),SNR_pdhg(1,end));
fprintf ('GRPDA:  %d/%f/%f\n',size(t_grpda,2)-1,t_grpda(1,end),SNR_grpda(1,end));
fprintf ('PCPDHG:  %d/%f/%f\n',size(t_pcpdhg,2)-1,t_pcpdhg(1,end),SNR_pcpdhg(1,end));
fprintf ('NPCPDHG:  %d/%f/%f\n',size(t_npcpdhg,2)-1,t_npcpdhg(1,end),SNR_npcpdhg(1,end));
fprintf ('DPHG:  %d/%f/%f\n',size(t_dphg,2)-1,t_dphg(1,end),SNR_dphg(1,end));












