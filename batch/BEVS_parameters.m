function par = BEVS_parameters()
    
    scaling=1e6; % % scaling parameter for numerical stability - don't change
    
    % cell growth and death
    mu=0.028; % 1/h
    k_deathT=8e-5; % uninfected cells death rate [1/h]
    k_deathI=0.0029; % infected cells death rate [1/h]
    tau_L=24; % infection age at which death kinetcs switches from k_deathT to  k_deathI [1/h]
    
    % binding and internalization
    tau_bind_decay=1.8;
    beta_decay=0.149;
    k_bind0=1.05e-8*60; % mL/cell/h 
    k_intern=0.01*60; % 1/h
    k_nucl=k_intern/2/scaling;
   
    % BV replication
    tau_repl_onset=6;
    tau_repl_end=18;
    k_repl=0.7318; %1/h

    % transcription parameters
    k_polh=607.89; % bp/h 
    k_p10=518.53; % bp/h
    tau_polh_onset=18; 
    tau_p10_onset=15; 
    tau_RNA_end=38; 

    % translation
    k_transl=2.94e9; 
    K_transl=2.74e4;

    % transgene amplification
    k_ampl=3.63e7; % #/cell/h 
    Kampl_rep=9.39e6; 

    % encapsidation
    k_pack=831.09; 
    Kencaps_coeff=5;    

    % degradation parameters
    kd_rep=0.03;% 1/h 
    kd_dna=0.01;%% 1/h 
    kd_gfp=3.2e-3; % 1/h 
    k_degrBV=7e-3; %1/h
    kd_rna=0.061; % 1/h

    % 3-Bac parameters
    tau_DIE1_onset=8;
    k_DIE1=k_polh/10;
    tau_rel_onset=18;
    k_rel=9.8;

    bp_vg=2700;
    bp_goi=700;

    %%
    par(1)=k_bind0;
    par(2)=mu;
    par(3)=k_deathT;
    par(4)=k_degrBV;
    par(5)=k_intern;
    par(6)=k_nucl;
    par(7)=k_repl;
    par(8)=tau_repl_end;
    par(9)=k_polh;
    par(10)=k_p10;
    par(11)=tau_polh_onset;
    par(12)=tau_p10_onset;
    par(13)=tau_RNA_end;
    par(14)=kd_rna;
    par(15)=k_transl;
    par(16)=K_transl;
    par(17)=k_pack;
    par(18)=Kampl_rep;
    par(19)=k_ampl;
    par(20)=kd_rep;
    par(21)=kd_dna;
    par(22)=kd_gfp;
    par(23)=tau_L;
    par(24)=k_deathI;
    par(25)=bp_vg;
    par(26)=bp_goi;
    par(27)=tau_bind_decay;
    par(28)=beta_decay;
    par(29)=tau_repl_onset;
    par(30)=Kencaps_coeff;
    par(31)=tau_DIE1_onset;
    par(32)=k_DIE1;
    par(33)=tau_rel_onset;
    par(34)=k_rel;
    par(35)=scaling;

end

