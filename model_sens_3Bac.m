function dxdt = model_sens_3Bac(t,x,par)

    k_bind0=par(1);
    mu=par(2);
    k_deathT=par(3);
    k_degrBV=par(4);
    k_intern=par(5);
    k_nucl=par(6);
    k_repl=par(7);
    tau_repl_end=par(8);
    k_polh=par(9);
    k_p10=par(10);
    tau_polh_onset=par(11);
    tau_p10_onset=par(12);
    tau_RNA_end=par(13);
    kd_rna=par(14);
    k_transl=par(15);
    K_transl=par(16);
    k_pack=par(17);
    Kgoi_rep=par(18);
    k_goi_repl=par(19);
    kd_rep=par(20);
    kd_dna=par(21);
    kd_gfp=par(22);
    tau_L=par(23);
    k_deathI=par(24);
    bp_goi=par(25);
    bp_goi_prod=par(26);
    tau_bind_decay=par(27);
    beta_decay=par(28);
    tau_repl_onset=par(29);
    K_rep_lim=par(30);

    % 3-Bac parameters
    tau_DIE1_onset=par(31);
    k_DIE1=par(32);
    tau_release_onset=par(33);
    release=par(34);
%     tau_release_end=par(35);
%     
%     BacN=3;
    Kgoi_dna=0.1; 
    bp_rep78=1863;
    bp_rep52=1191;
    bp_cap=1557; 
    scaling=1e6; 

    dxdt=zeros(length(x),1);
    kd_cap=0;
    tot_cells_1stw=sum(x(5:11));
    
    % binding first wave
    if t>tau_bind_decay
        bind_decay1=exp(-beta_decay*(t-tau_bind_decay));
    else
        bind_decay1=1;
    end
    k_bind1=k_bind0*bind_decay1;
    k_bind01=k_bind0*1;

    if t>15
        k_bind1=0;
        k_bind01=0;
    end
    
    %% death rates calculation
    % cells death rate transition 
    k_death1=zeros(1,7);
    for i=1:7
        if x(4+i)>1e-6
            logDNA=log((x(15+(i-1)*22)+x(16+(i-1)*22)+...
                    x(17+(i-1)*22))/x(4+i)*scaling+1e-16);
            if logDNA<1
                logDNA=1;
            end
            k_death1(i)=k_deathT*(0.5+0.5*tanh((tau_L-t)/0.3))+...
                k_deathI*logDNA.*(0.5+0.5*tanh((t-tau_L)/0.3));
        else
            k_death1(i)=k_deathT;
        end
    end
    
    %% release calculation (1st wave)
    release_rep=0;
    release_cap=0;
    release_goi=0;
    if t>tau_release_onset  %&& t<tau_release_end
        for i=1:7
            if x(4+i)>0
                total_dna=x(15+(i-1)*22)+x(16+(i-1)*22)+x(17+(i-1)*22);
                if total_dna>0
                    release_rep=release_rep+x(4+i)*release*x(15+(i-1)*22)/total_dna;
                    release_cap=release_cap+x(4+i)*release*x(16+(i-1)*22)/total_dna;
                    release_goi=release_goi+x(4+i)*release*x(17+(i-1)*22)/total_dna;
                end
            end
        end
    end

    %% Balances
    %% virions balances
    dxdt(1)=mu*x(1)-k_deathT*x(1)-k_bind01*(x(1)*x(2))-...
        k_bind01*(x(1)*x(3))-k_bind01*(x(1)*x(4));

    dxdt(2)=-k_bind01*x(2)*x(1)-k_bind1*x(2)*tot_cells_1stw-...
        k_degrBV*x(2)+release_rep; % repBV (#/mL) 
    dxdt(3)=-k_bind01*x(3)*x(1)-k_bind1*x(3)*tot_cells_1stw-...
        -k_degrBV*x(3)+release_cap; % capBV (#/mL) 
    dxdt(4)=-k_bind01*x(4)*x(1)-k_bind1*x(4)*tot_cells_1stw-...
        -k_degrBV*x(4)+release_goi; % goiBV (#/mL) 
    
    %% cells first wave balances
    dxdt(5)=k_bind01*x(1)*x(2)-k_bind1*(x(5)*(x(3)+x(4)))-k_death1(1)*x(5); % infected by repBV (#/mL)
    dxdt(6)=k_bind01*x(1)*x(3)-k_bind1*(x(6)*(x(2)+x(4)))-k_death1(2)*x(6); % infected by capBV (#/mL)
    dxdt(7)=k_bind01*x(1)*x(4)-k_bind1*(x(7)*(x(2)+x(3)))-k_death1(3)*x(7); % infected by goiBV (#/mL)
   
    dxdt(8)=k_bind1*(x(5)*x(3)+x(6)*x(2)-x(8)*x(4))-k_death1(4)*x(8); % infected by repBV and capBV (#/mL)
    dxdt(9)=k_bind1*(x(5)*x(4)+x(7)*x(2)-x(9)*x(3))-k_death1(5)*x(9); % infected by repBV and goiBV (#/mL)
    dxdt(10)=k_bind1*(x(6)*x(4)+x(7)*x(3)-x(10)*x(2))-k_death1(6)*x(10); % infected by capBV and goiBV (#/mL)
    
    dxdt(11)=k_bind1*(x(8)*x(4)+x(9)*x(3)+x(10)*x(2))-k_death1(7)*x(11); % triple infected
    
    %% virus bound per infected cell (#/mL)
    % infected by one BV
    dxdt(12)=k_bind01*x(2)*x(1)+k_bind1*x(2)*x(5)-k_bind1*x(12)*(x(3)+x(4))-...
        k_intern*x(12)-k_death1(1)*x(12);
    dxdt(13)=0;
    dxdt(14)=0;

    dxdt(34)=0;
    dxdt(35)=k_bind01*x(3)*x(1)+k_bind1*x(3)*x(6)-k_bind1*x(35)*(x(2)+x(4))-...
        k_intern*x(35)-k_death1(2)*x(35);
    dxdt(36)=0;

    dxdt(56)=0;
    dxdt(57)=0;
    dxdt(58)=k_bind01*x(4)*x(1)+k_bind1*x(4)*x(7)-k_bind1*x(58)*(x(2)+x(3))-...
        k_intern*x(58)-k_death1(3)*x(58);

    % infected by two BV 
    dxdt(78)=k_bind1*x(2)*(x(6)+x(8))+k_bind1*x(12)*x(3)-... % rep_cap - x(8)
        k_bind1*x(78)*x(4)-k_intern*x(78)-k_death1(4)*x(78);
    dxdt(79)=k_bind1*x(3)*(x(5)+x(8))+k_bind1*x(35)*x(2)-...
        k_bind1*x(79)*x(4)-k_intern*x(79)-k_death1(4)*x(79);
    dxdt(80)=0;

    dxdt(100)=k_bind1*x(2)*(x(7)+x(9))+k_bind1*x(12)*x(4)-... % rep_goi - x(9)
        k_bind1*x(100)*x(3)-k_intern*x(100)-k_death1(5)*x(100);
    dxdt(101)=0;
    dxdt(102)=k_bind1*x(4)*(x(5)+x(9))+k_bind1*x(58)*x(2)-...
        k_bind1*x(102)*x(3)-k_intern*x(102)-k_death1(5)*x(102);

    dxdt(122)=0;                                            % cap_goi - x(10)
    dxdt(123)=k_bind1*x(3)*(x(7)+x(10))+k_bind1*x(35)*x(4)-... 
        k_bind1*x(123)*x(2)-k_intern*x(123)-k_death1(6)*x(123);
    dxdt(124)=k_bind1*x(4)*(x(6)+x(10))+k_bind1*x(58)*x(3)-...
        k_bind1*x(124)*x(2)-k_intern*x(124)-k_death1(6)*x(124);
    
    % triple infected - x(11)
    dxdt(144)=k_bind1*x(2)*(x(10)+x(11))+...
        k_bind1*x(78)*x(4)+k_bind1*x(100)*x(3)-k_intern*x(144)-k_death1(7)*x(144);
    
    dxdt(145)=k_bind1*x(3)*(x(9)+x(11))+...
        k_bind1*x(79)*x(4)+k_bind1*x(123)*x(2)-k_intern*x(145)-k_death1(7)*x(145);
    
    dxdt(146)=k_bind1*x(4)*(x(8)+x(11))+...
        k_bind1*x(102)*x(3)+k_bind1*x(124)*x(2)-k_intern*x(146)-k_death1(7)*x(146);
    
    %% virus in nucleus: transport

    % infected by one BV
    dxdt(15)=-k_bind1*x(15)*(x(3)+x(4))+k_nucl*x(12)-k_death1(1)*x(15); 
    dxdt(16)=0;
    dxdt(17)=0;

    dxdt(37)=0;
    dxdt(38)=-k_bind1*x(38)*(x(2)+x(4))+k_nucl*x(35)-k_death1(2)*x(38);
    dxdt(39)=0;

    dxdt(59)=0;
    dxdt(60)=0;
    dxdt(61)=-k_bind1*x(61)*(x(2)+x(3))+k_nucl*x(58)-k_death1(3)*x(61);

    % infected by two BV 
    dxdt(81)=k_bind1*x(15)*x(3)-... % rep_cap - x(8)
        k_bind1*x(81)*x(4)+k_nucl*x(78)-k_death1(4)*x(81);
    dxdt(82)=k_bind1*x(38)*x(2)-...
        k_bind1*x(82)*x(4)+k_nucl*x(79)-k_death1(4)*x(82);
    dxdt(83)=0;

    dxdt(103)=k_bind1*x(15)*x(4)-... % rep_goi - x(9)
        k_bind1*x(103)*x(3)+k_nucl*x(100)-k_death1(5)*x(103);
    dxdt(104)=0;
    dxdt(105)=k_bind1*x(61)*x(2)-...
        k_bind1*x(105)*x(3)+k_nucl*x(102)-k_death1(5)*x(105);

    dxdt(125)=0;
    dxdt(126)=k_bind1*x(38)*x(4)-... % cap_goi - x(10)
        k_bind1*x(126)*x(2)+k_nucl*x(123)-k_death1(6)*x(126);
    dxdt(127)=k_bind1*x(61)*x(3)-...
        k_bind1*x(127)*x(2)+k_nucl*x(124)-k_death1(6)*x(127);

    % triple infected
    dxdt(147)=k_bind1*x(81)*x(4)+k_bind1*x(103)*x(3)+k_nucl*x(144)-...
        k_death1(7)*x(147);
    dxdt(148)=k_bind1*x(82)*x(4)+k_bind1*x(126)*x(2)+k_nucl*x(145)-...
        k_death1(7)*x(148);
    dxdt(149)=k_bind1*x(105)*x(3)+k_bind1*x(127)*x(2)+k_nucl*x(146)-...
        k_death1(7)*x(149);
    %%%%%%%%

    %% cell content
    % replication, transcription and translation activation functions
    f_repl=(0.5+0.5*tanh((t-tau_repl_onset)/0.3)).*(0.5+0.5*tanh((tau_repl_end-t)/1)); % smooth transition

    f_DIE1=(0.5+0.5*tanh((t-tau_DIE1_onset)/0.3)).*(0.5+0.5*tanh((tau_RNA_end-t)/1.1)); % smooth transition
    f_polh=(0.5+0.5*tanh((t-tau_polh_onset)/0.3)).*(0.5+0.5*tanh((tau_RNA_end-t)/1.1)); % smooth transition
    f_p10=(0.5+0.5*tanh((t-tau_p10_onset)/0.3)).*(0.5+0.5*tanh((tau_RNA_end-t)/1.1)); % smooth transition

    for i =1:7
        if x(i+4)>1e-16
            % viral DNA replication
            dxdt(15+(i-1)*22)=dxdt(15+(i-1)*22)+k_repl*f_repl*x(15+(i-1)*22); % nuclear rep DNA
            dxdt(16+(i-1)*22)=dxdt(16+(i-1)*22)+k_repl*f_repl*x(16+(i-1)*22);  % nuclear cap DNA
            dxdt(17+(i-1)*22)=dxdt(17+(i-1)*22)+k_repl*f_repl*x(17+(i-1)*22); % nuclear GOI DNA
            
            % transcription (3-Bac system) - #/mL (scaled)
%             if BacN==3
            dxdt(18+(i-1)*22)=k_DIE1/bp_rep78*f_DIE1*x(15+(i-1)*22)-...
                kd_rna*(x(18+(i-1)*22))-k_death1(i)*x(18+(i-1)*22); % rep78 mRNA
            dxdt(19+(i-1)*22)=k_polh/bp_rep52*f_polh*x(15+(i-1)*22)-...
                kd_rna*(x(19+(i-1)*22))-k_death1(i)*x(19+(i-1)*22); % rep52 mRNA
            dxdt(20+(i-1)*22)=k_polh/bp_cap*f_polh*x(16+(i-1)*22)-...
                kd_rna*x(20+(i-1)*22)-k_death1(i)*x(20+(i-1)*22); % VP mRNA
            dxdt(21+(i-1)*22)=k_p10/bp_goi_prod*f_p10*(x(17+(i-1)*22)+x(22+(i-1)*22))-...
                kd_rna*(x(21+(i-1)*22))-k_death1(i)*x(21+(i-1)*22); % transgene mRNA
%             else
%                 % transcription (2-Bac system) - #/mL (scaled) 
%                 dxdt(19+(i-1)*22)=k_polh/bp_rep52*f_polh*x(15+(i-1)*22)-...
%                     kd_rna*(x(19+(i-1)*22))-k_death1(i)*x(19+(i-1)*22); % rep52 mRNA
%                 dxdt(20+(i-1)*22)=k_p10/bp_cap*f_p10*x(15+(i-1)*22)-...
%                     kd_rna*x(20+(i-1)*22)-k_death1(i)*x(20+(i-1)*22); % VP mRNA
%                 dxdt(21+(i-1)*22)=k_p10/bp_goi_prod*f_p10*(x(17+(i-1)*22)+x(22+(i-1)*22))-...
%                     kd_rna*(x(21+(i-1)*22))-k_death1(i)*x(21+(i-1)*22); % transgene mRNA
%             end

            % transgene amplification - #/mL (scaled)
            dxdt(22+(i-1)*22)=k_goi_repl*2700/bp_goi/scaling*x(i+4)*...
                (x(17+(i-1)*22)+x(22+(i-1)*22))/...
                (Kgoi_dna/scaling*x(i+4)+x(17+(i-1)*22)+x(22+(i-1)*22))*...
                x(23+(i-1)*22)/...
                (Kgoi_rep/scaling*x(i+4)+x(23+(i-1)*22))-kd_dna*x(22+(i-1)*22)-...
                k_death1(i)*x(22+(i-1)*22);
            
            % translation - #/mL (scaled)
            r(1)=k_transl/scaling*x(i+4)*(x(18+(i-1)*22))/(K_transl/...
                scaling*x(i+4)+x(18+(i-1)*22));
            r(2)=k_transl/scaling*x(i+4)*(x(19+(i-1)*22))/(K_transl/...
                scaling*x(i+4)+x(19+(i-1)*22));
            r(3)=k_transl/scaling*x(i+4)*(x(20+(i-1)*22))/(K_transl/...
                scaling*x(i+4)+x(20+(i-1)*22));
            r(4)=k_transl/scaling*x(i+4)*(x(21+(i-1)*22))/(K_transl/...
                scaling*x(i+4)+x(21+(i-1)*22));

%             if BacN==2
%                 r(2)=r(2)/2;
%                 r(1)=r(2);
%                 bp_rep78=bp_rep52;
%             end
            
            if sum(r)>k_transl*x(i+4)/scaling
                r=r./sum(r)*k_transl*x(i+4)/scaling;
            end

            dxdt(23+(i-1)*22)=r(1)/bp_rep78-kd_rep*(x(23+(i-1)*22))-k_death1(i)*x(23+(i-1)*22); % rep78
            dxdt(24+(i-1)*22)=r(2)/bp_rep52-kd_rep*(x(24+(i-1)*22))-k_death1(i)*x(24+(i-1)*22); % rep52
            dxdt(25+(i-1)*22)=r(3)/bp_cap/60-k_death1(i)*x(25+(i-1)*22)-kd_cap*(x(25+(i-1)*22)); % capsids
            dxdt(26+(i-1)*22)=r(4)/bp_goi_prod-kd_gfp*(x(26+(i-1)*22))-k_death1(i)*x(26+(i-1)*22); % GFP
            
            % packaging
            % 1st order packaging with rep limitation + limiting species
            if x(22+(i-1)*22)>x(25+(i-1)*22)
                pack_lim=x(25+(i-1)*22);
            else
                pack_lim=x(22+(i-1)*22);
            end
            K_rep=pack_lim*K_rep_lim+1;
            r_pack=k_pack/bp_goi*pack_lim*...
                x(24+(i-1)*22)/(x(24+(i-1)*22)+K_rep/scaling*x(i+4));

            dxdt(25+(i-1)*22)=dxdt(25+(i-1)*22)-r_pack;
            dxdt(22+(i-1)*22)=dxdt(22+(i-1)*22)-r_pack;
            dxdt(27+(i-1)*22)=r_pack-kd_cap*(x(27+(i-1)*22)); 

            dxdt(28+(i-1)*22) = k_death1(i)*x(4+i); % number of dead cells
            dxdt(29+(i-1)*22) = k_death1(i)*x(22+(i-1)*22)-...
                kd_dna*x(29+(i-1)*22); % number of GOI copies in dead cells
            dxdt(30+(i-1)*22) = k_death1(i)*x(23+(i-1)*22)-kd_rep*(x(30+(i-1)*22)); % rep78 protein conc.
            dxdt(31+(i-1)*22) = k_death1(i)*x(24+(i-1)*22)-kd_rep*(x(31+(i-1)*22)); % rep52 protein conc.
            dxdt(32+(i-1)*22) = k_death1(i)*x(25+(i-1)*22)-kd_cap*(x(32+(i-1)*22)); % empty capsids conc.
        end
    end
end
