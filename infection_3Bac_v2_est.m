function [t,x] = infection_3Bac_v2_est(par,MOI, k_deathI, tau_L, C0, ...
    bp_goi, bp_goi_prod, BacN,S,t_end)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
% x(1) = target cell (#/mL)
%
% virions
% x(2) = BV rep (#/mL) 
% x(3) = BV cap (#/mL) 
% x(4) = BV GOI (#/mL)
%
% cells infected by one virus
% x(5) = rep
% x(6) = cap
% x(7) = goi
%
% cells infected by two viruses
% x(8) = rep_cap
% x(9) = rep_GOI
% x(10) = cap_goi
% 
% cells infected by three viruses
% x(11) = triple-infected
%
% content of infected cells (for i=1:7)
% x(12+(i-1)*22) = bound repBV
% x(13+(i-1)*22) = bound capBV
% x(14+(i-1)*22) = bound GOIBV
% x(15+(i-1)*22) = nuclear rep DNA
% x(16+(i-1)*22) = nuclear cap DNA
% x(17+(i-1)*22) = nuclear GOI DNA
% x(18+(i-1)*22) = rep78 mRNA
% x(19+(i-1)*22) = rep52 mRNA
% x(20+(i-1)*22) = VP mRNA
% x(21+(i-1)*22) = GOI mRNA
% x(22+(i-1)*22) = GOI number of copies
% x(23+(i-1)*22) = rep78 protein conc.
% x(24+(i-1)*22) = rep52 protein conc.
% x(25+(i-1)*22) = empty capsids conc.
% x(26+(i-1)*20)) = GOI protein
% x(27+(i-1)*22) = filled capsids conc.
% 
% content of dead cells
% x(28+(i-1)*22)) = number of dead triple-infected cells
% x(29+(i-1)*22)) = number of GOI copies in dead cells
% x(30+(i-1)*22)) = rep78 protein conc.
% x(31+(i-1)*22)) = rep52 protein conc.
% x(32+(i-1)*22)) = empty capsids conc.
% x(33+(i-1)*22)) = filled capsids conc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 10
    t_end=72;
end

MOIrep=MOI(2);

if BacN==3
    MOIcap=MOI(3);
elseif BacN==2
    MOIcap=0;
end
MOIgoi=MOI(1);

Vrep0=C0*MOIrep; % #/mL
Vcap0=C0*MOIcap; % #/mL
Vgoi0=C0*MOIgoi; % #/mL

x0=[C0 Vrep0 Vcap0 Vgoi0 zeros(1,360)];

options=odeset('RelTol',1e-4,'JPattern',S);

[k_bind0,mu, k_deathT, k_degrBV, k_intern, k_nucl, k_repl,...
        ~, tau_repl_end, k_polh, k_p10, k_DIE1, tau_DIE1_start,...
        tau_polh_start, tau_p10_start,tau_RNA_end, kd_rna, k_transl, ...
        K_transl, k_pack, bp_rep78, bp_rep52, bp_cap, cap_eff, rep_eff, Kgoi_dna,...
        Kgoi_rep, k_goi_repl, ~, ~, ~, kd_rep,...
        kd_cap, kd_dna, kd_gfp, scaling] = BEVS_parameters_est(par);

[t,x]=ode45(@inf_model,0:.1:t_end,x0,options,k_bind0,mu, k_deathT, k_degrBV, k_intern, k_nucl, k_repl,...
        tau_repl_end, k_polh, k_p10, k_DIE1, tau_DIE1_start,...
        tau_polh_start, tau_p10_start,tau_RNA_end, kd_rna, k_transl, ...
        K_transl, k_pack, bp_rep78, bp_rep52, bp_cap, cap_eff, rep_eff, Kgoi_dna,...
        Kgoi_rep, k_goi_repl, kd_rep,...
        kd_cap, kd_dna, kd_gfp, scaling, BacN, tau_L, k_deathI, bp_goi, bp_goi_prod);

end