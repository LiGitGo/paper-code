%% load data
% load('data/yellow_flower_beta10_alpha1_tau1_500s_pdhg_1.mat');
% load('data/blue_flower_300_beta10_alpha1_tau1_500s_pdhg_1.mat');

%% 

%Tol = 1e-7
index_PDA_e7 = find(-((f_pdhg-optval)/optval)< 1e-7, 1);
index_PC_e7 = find(-((f_pcpdhg-optval)/optval)< 1e-7, 1);
index_NPC_e7 = find(-((f_npcpdhg-optval)/optval)< 1e-7, 1);
index_GNPC_e7 = find(-((f_dphg(11:end)-optval)/optval)< 1e-7, 1) + 10;
index_GRPDA_e7 = find(-((f_grpda-optval)/optval)< 1e-7, 1);

time_PDA_e7 = t_pdhg(index_PDA_e7);
time_PC_e7 = t_pcpdhg(index_PC_e7);
time_NPC_e7 = t_npcpdhg(index_NPC_e7);
time_GNPC_e7 = t_dphg(index_GNPC_e7);
time_GRPDA_e7 = t_grpda(index_GRPDA_e7);

fprintf ('Image=%s,  Tol=%e\n',name,1e-7);
fprintf ('PDHG:  %d/%f\n',index_PDA_e7,time_PDA_e7);
fprintf ('GRPDA:  %d/%f\n',index_GRPDA_e7,time_GRPDA_e7);
fprintf ('PCPDHG:  %d/%f\n',index_PC_e7,time_PC_e7);
fprintf ('NPCPDHG:  %d/%f\n',index_NPC_e7,time_NPC_e7);
fprintf ('DPHG:  %d/%f\n',index_GNPC_e7,time_GNPC_e7);
fprintf ('\n');

%Tol = 1e-8
index_PDA_e8 = find(-((f_pdhg-optval)/optval)< 1e-8, 1);
index_PC_e8 = find(-((f_pcpdhg-optval)/optval)< 1e-8, 1);
index_NPC_e8 = find(-((f_npcpdhg-optval)/optval)< 1e-8, 1);
index_GNPC_e8 = find(-((f_dphg(11:end)-optval)/optval)< 1e-8, 1) + 10;
index_GRPDA_e8 = find(-((f_grpda-optval)/optval)< 1e-8, 1);

time_PDA_e8 = t_pdhg(index_PDA_e8);
time_PC_e8 = t_pcpdhg(index_PC_e8);
time_NPC_e8 = t_npcpdhg(index_NPC_e8);
time_GNPC_e8 = t_dphg(index_GNPC_e8);
time_GRPDA_e8 = t_grpda(index_GRPDA_e8);

fprintf ('Image=%s,  Tol=%e\n',name,1e-8);
fprintf ('PDHG:  %d/%f\n',index_PDA_e8,time_PDA_e8);
fprintf ('GRPDA:  %d/%f\n',index_GRPDA_e8,time_GRPDA_e8);
fprintf ('PCPDHG:  %d/%f\n',index_PC_e8,time_PC_e8);
fprintf ('NPCPDHG:  %d/%f\n',index_NPC_e8,time_NPC_e8);
fprintf ('DPHG:  %d/%f\n',index_GNPC_e8,time_GNPC_e8);
fprintf ('\n');

%Tol = 1e-9
index_PDA_e9 = find(-((f_pdhg-optval)/optval)< 1e-9, 1);
index_PC_e9 = find(-((f_pcpdhg-optval)/optval)< 1e-9, 1);
index_NPC_e9 = find(-((f_npcpdhg-optval)/optval)< 1e-9, 1);
index_GNPC_e9 = find(-((f_dphg(11:end)-optval)/optval)< 1e-9, 1) + 10;
index_GRPDA_e9 = find(-((f_grpda-optval)/optval)< 1e-9, 1);

time_PDA_e9 = t_pdhg(index_PDA_e9);
time_PC_e9 = t_pcpdhg(index_PC_e9);
time_NPC_e9 = t_npcpdhg(index_NPC_e9);
time_GNPC_e9 = t_dphg(index_GNPC_e9);
time_GRPDA_e9 = t_grpda(index_GRPDA_e9);

fprintf ('Image=%s,  Tol=%e\n',name,1e-9);
fprintf ('PDHG:  %d/%f\n',index_PDA_e9,time_PDA_e9);
fprintf ('GRPDA:  %d/%f\n',index_GRPDA_e9,time_GRPDA_e9);
fprintf ('PCPDHG:  %d/%f\n',index_PC_e9,time_PC_e9);
fprintf ('NPCPDHG:  %d/%f\n',index_NPC_e9,time_NPC_e9);
fprintf ('DPHG:  %d/%f\n',index_GNPC_e9,time_GNPC_e9);

