function [k_bind0,mu, k_deathT, k_degrBV, k_intern, k_nucl, k_repl,...
    tau_VL, tau_repl_end, k_polh, k_p10, k_DIE1, tau_DIE1_start,...
    tau_polh_start, tau_p10_start,tau_RNA_end, kd_rna, k_transl, ...
    K_transl, k_pack, bp_rep78, bp_rep52, bp_cap, cap_eff, rep_eff, Kgoi_dna,...
    Kgoi_rep, k_goi_repl, Kpack_rep, Kpack_goi, Kpack_cap, kd_rep,...
    kd_cap, kd_dna, kd_gfp, scaling] = BEVS_parameters_est(par)

    k_transl=10^par(1);
    K_transl=10^par(2);
    
    Kgoi_rep=10.^par(3); %11e5;
    k_goi_repl=10.^par(4); %12500; % #/cell/h (from Tam2) 
    k_pack=10.^par(5); %103*ATP_conc/(31+ATP_conc)*10; % bp/h - max velocity of vg encapsidation

    rep_eff=1;

    scaling=1e6;

    Kgoi_dna=0.1;
    k_bind0=1.05e-8*60; % mL/cell/h Nielsen 2003
    mu=0.028; % 1/h
    k_deathT=8e-5; % 1/h
    k_degrBV=7e-3; %1/h
    k_intern=2e-3*60; % 1/h
    k_nucl=k_intern/2/scaling;
    
    k_repl=0.7303; %1/h
    tau_VL=18; % h
    tau_repl_end=tau_VL;
    
    % transcription parameters
    k_polh=607.89; % bp/h - cap
    k_p10=518.53; % bp/h 
    
    tau_DIE1_start=8;
    tau_polh_start=18; %20;
    tau_p10_start=15; 
    tau_RNA_end=38; %x(4);
    
    kd_rna=0.061; %0.05; %x(2); % degradation rate
    
    bp_rep78=1863;
    bp_rep52=1191;
    bp_cap=1557; % single transcript (vp1)

    % uncertain parameters
    k_DIE1=k_polh/10; % bp/h - rep78
    cap_eff=1;

    Kpack_rep=1000;%6e6; %150000.11*1e6; % #/cell (for rep from Tam2) 
    Kpack_goi=1000;%00;%6.11e5; % #/cell (I'm assuming it)
    Kpack_cap=1000;%0.11e3; % #/cel (I'm assuming it)
    
    kd_rep=0.03;%7.45e-2; % 1/h (from Tam2)
    kd_cap=0;%6.28e-2; % 1/h (from Tam2)
    kd_dna=0.01;%6.28e-2; % 1/h (from Tam2)
    kd_gfp=3.2e-3; % 1/h (I'm assuming it)

        p=[9.4662    4.4358    6.9768    4.1389    2.9165];

    scaling=1e6;
    
    % cell growth and death
    mu=0.028; % 1/h
    k_deathT=8e-5; % 1/h
    k_deathI=0.0029;
    tau_L=24;
    
    % binding and internalization
    tau_bind_decay=1.8;
    beta_decay=0.149;
    k_bind0=1.05e-8*60; % mL/cell/h Nielsen 2003
    k_intern=0.01*60; % 1/h
    k_nucl=k_intern/2/scaling;
   
    % BV replication
    tau_repl_onset=6;
    tau_repl_end=18;
    k_repl=0.7303; %1/h

    % transcription parameters
    k_polh=607.89; % bp/h - cap
    k_p10=518.53; % bp/h
    tau_polh_onset=18; %20;
    tau_p10_onset=15; 
    tau_RNA_end=38; %x(4);
    
    % translation
    k_transl=10^p(1); %10^par(1);
    K_transl=10^p(2);

    % transgene amplification
    Kgoi_rep=10.^p(3); %11e5;
    k_goi_repl=10.^p(4); %12500; % #/cell/h (from Tam2) 

    % encapsidation
    K_rep_lim=5;    
    k_pack=10.^p(5); %103*ATP_conc/(31+ATP_conc)*10; % bp/h - max velocity of vg encapsidation

    % degradation parameters
    kd_rep=0.03;%7.45e-2; % 1/h (from Tam2)
    kd_dna=0.01;%6.28e-2; % 1/h (from Tam2)
    kd_gfp=3.2e-3; % 1/h (I'm assuming it)
    k_degrBV=7e-3; %1/h
    kd_rna=0.061; % degradation rate

    % 3-Bac parameters
    tau_DIE1_onset=8;
    k_DIE1=k_polh/10;
    tau_release_onset=18;
    release=9.8;
    tau_release_end=45;

    bp_goi=2700;
    bp_goi_prod=700;

end

 