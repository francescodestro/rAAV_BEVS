function dxdt = BEVS_model(t,x,BacN,par)
    %% parameters
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
    tau_polh_on=par(11);
    tau_p10_on=par(12);
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
    tau_repl_on=par(29);
    K_encaps_coeff=par(30);

    % 3-Bac parameters
    tau_DIE1_on=par(31);
    k_DIE1=par(32);
    tau_rel_on=par(33);
    k_rel=par(34);
    scaling=par(35); 

    Kgoi_dna=0.1; 
    bp_rep78=1863;
    bp_rep52=1191;
    
    if BacN==2
        bp_avg=0.5*(bp_rep78+bp_rep52);
        bp78_transl=bp_avg;
        bp52_transl=bp_avg;
    else
        bp78_transl=bp_rep78;
        bp52_transl=bp_rep52;
    end
    bp_cap=1557; 


    %% simulation
    dxdt=zeros(length(x),1);
    x(x<0)=0; % for stability

    tot_cells_1stw=sum(x(5:11));
    tot_cells_2ndw=sum(x(205:211));
    
    % binding first wave
    if t>1.8
        bind_decay1=exp(-beta_decay*(t-tau_bind_decay));
    else
        bind_decay1=1;
    end
    k_bind1=k_bind0*bind_decay1;

    % binding second wave
    if t>tau_rel_on
        if t>tau_rel_on+tau_bind_decay
            bind_decay2=exp(-beta_decay*(t-tau_rel_on-tau_bind_decay));
        else
            bind_decay2=1;
        end
        k_bind1=0;
        k_bind01=0;
        k_bind02=k_bind0;
        k_bind2=k_bind02*bind_decay2;

    else
        k_bind01=k_bind0;
        k_bind2=0;
        k_bind02=0;
    end

    %% death rates calculation
    % cells death rate transition 
    k_death1=zeros(1,7);
    k_death2=zeros(1,7);
    for i=1:7
        if x(4+i)>1e-6
            k_death1(i)=k_deathT*(0.5+0.5*tanh((tau_L-t)/0.3))+...
                (k_deathI(1)*max(1,log((x(15+(i-1)*22)+x(16+(i-1)*22)+...
                    x(17+(i-1)*22))/x(4+i)*scaling+1e-16))+...
                0*max(0,log(x(23+(i-1)*22)/x(4+i)*scaling+1e-16))+...
                0*max(0,log((x(25+(i-1)*22)+x(27+(i-1)*22))/x(4+i)*scaling+1e-16))).*...
                (0.5+0.5*tanh((t-tau_L)/0.3));
        else
            k_death1(i)=k_deathT;
        end
        if x(204+i)>1e-6
            k_death2(i)=k_deathT*(0.5+0.5*tanh((tau_L+tau_rel_on-t)/0.3))+...
                (k_deathI(1)*max(1,log((x(215+(i-1)*22)+x(216+(i-1)*22)+...
                    x(217+(i-1)*22)/x(204+i)*scaling+1e-16)))+...
                0*max(0,log(x(223+(i-1)*22)/x(204+i)*scaling+1e-16))+...
                0*max(0,log((x(225+(i-1)*22)+x(227+(i-1)*22))/x(204+i)*scaling+1e-16))).*...
                (0.5+0.5*tanh((t-tau_L-18)/0.3));
        else
            k_death2(i)=k_deathT;
        end
    end

    %% release calculation (1st wave)
    release_rep=0;
    release_cap=0;
    release_goi=0;
    if t>tau_rel_on 
        for i=1:7
            if x(4+i)>0
                total_dna=x(15+(i-1)*22)+x(16+(i-1)*22)+x(17+(i-1)*22);
                if total_dna>0
                    release_rep=release_rep+x(4+i)*k_rel*x(15+(i-1)*22)/total_dna;
                    release_cap=release_cap+x(4+i)*k_rel*x(16+(i-1)*22)/total_dna;
                    release_goi=release_goi+x(4+i)*k_rel*x(17+(i-1)*22)/total_dna;
                end
            end
            if t>tau_rel_on 
                if x(204+i)>0
                    total_dna=x(215+(i-1)*22)+x(216+(i-1)*22)+x(217+(i-1)*22);
                    if total_dna>0
                        release_rep=release_rep+x(204+i)*k_rel*x(215+(i-1)*22)/total_dna;
                        release_cap=release_cap+x(204+i)*k_rel*x(216+(i-1)*22)/total_dna;
                        release_goi=release_goi+x(204+i)*k_rel*x(217+(i-1)*22)/total_dna;
                    end
                end
            end
        end
    end

    %% Balances
    %% virions balances
    dxdt(1)=mu*x(1)-k_deathT*x(1)-k_bind0*(x(1)*x(2))-...
        k_bind0*(x(1)*x(3))-k_bind0*(x(1)*x(4));

    dxdt(2)=-k_bind0*x(2)*x(1)-k_bind1*x(2)*tot_cells_1stw-...
        k_bind2*x(2)*tot_cells_2ndw-k_degrBV*x(2)+release_rep; % repcapBV (TwoBac) or repBV (ThreeBac) [#/mL] 
    dxdt(3)=-k_bind0*x(3)*x(1)-k_bind1*x(3)*tot_cells_1stw-...
        k_bind2*x(3)*tot_cells_2ndw-k_degrBV*x(3)+release_cap; % capBV (ThreeBac) [#/mL] 
    dxdt(4)=-k_bind0*x(4)*x(1)-k_bind1*x(4)*tot_cells_1stw-...
        k_bind2*x(4)*tot_cells_2ndw-k_degrBV*x(4)+release_goi; % goiBV (TwoBac and ThreeBac) [#/mL] 
    
    %% cells first wave balances
    dxdt(5)=k_bind01*x(1)*x(2)-k_bind1*(x(5)*(x(3)+x(4)))-k_death1(1)*x(5); % infected by repcapBV (TwoBac) or repBV (ThreeBac) [#/mL]
    dxdt(6)=k_bind01*x(1)*x(3)-k_bind1*(x(6)*(x(2)+x(4)))-k_death1(2)*x(6); % infected by capBV (ThreeBac) [#/mL]
    dxdt(7)=k_bind01*x(1)*x(4)-k_bind1*(x(7)*(x(2)+x(3)))-k_death1(3)*x(7); % infected by goiBV [#/mL]
   
    dxdt(8)=k_bind1*(x(5)*x(3)+x(6)*x(2)-x(8)*x(4))-k_death1(4)*x(8); % infected by repBV and capBV (ThreeBac) [#/mL]
    dxdt(9)=k_bind1*(x(5)*x(4)+x(7)*x(2)-x(9)*x(3))-k_death1(5)*x(9); % infected by repcapBV and goiBV (TwoBac) or by repBV and goiBV (ThreeBac) (#/mL)
    dxdt(10)=k_bind1*(x(6)*x(4)+x(7)*x(3)-x(10)*x(2))-k_death1(6)*x(10); % infected by capBV and goiBV (ThreeBac) [#/mL]
    
    dxdt(11)=k_bind1*(x(8)*x(4)+x(9)*x(3)+x(10)*x(2))-k_death1(7)*x(11); % infected by repBV, capBV, goiBV (ThreeBac) [#/mL]
    
    %% virus bound to infected cells (#/mL)

    % in cells infected by one BV:
        % by repcapBV (TwoBac) or repBV (ThreeBac)
        dxdt(12)=k_bind01*x(2)*x(1)+k_bind1*x(2)*x(5)-k_bind1*x(12)*(x(3)+x(4))-...
            k_intern*x(12)-k_death1(1)*x(12); % receptor-bound repcapBV (TwoBac) or repBV (ThreeBac) [#/mL]
        dxdt(13)=0; % receptor-bound capBV (ThreeBac) [#/mL]
        dxdt(14)=0; % receptor-bound goiBV [#/mL]
        
        % by capBV (ThreeBac) 
        dxdt(34)=0; % receptor-bound repcapBV (TwoBac) or repBV (ThreeBac) [#/mL]
        dxdt(35)=k_bind01*x(3)*x(1)+k_bind1*x(3)*x(6)-k_bind1*x(35)*(x(2)+x(4))-...
            k_intern*x(35)-k_death1(2)*x(35); % receptor-bound capBV (ThreeBac) [#/mL]
        dxdt(36)=0; %  receptor-bound goiBV [#/mL]
        
        % by goiBV:
        dxdt(56)=0; % receptor-bound repcapBV (TwoBac) or repBV (ThreeBac) [#/mL]
        dxdt(57)=0; % receptor-bound capBV (ThreeBac) [#/mL]
        dxdt(58)=k_bind01*x(4)*x(1)+k_bind1*x(4)*x(7)-k_bind1*x(58)*(x(2)+x(3))-...
            k_intern*x(58)-k_death1(3)*x(58); % receptor-bound goiBV [#/mL]

    % in cells infected by two BV:
        % infected by repBV and capBV (ThreeBac) [#/mL]
        dxdt(78)=k_bind1*x(2)*(x(6)+x(8))+k_bind1*x(12)*x(3)-... % receptor-bound repcapBV (TwoBac) or repBV (ThreeBac) [#/mL]
            k_bind1*x(78)*x(4)-k_intern*x(78)-k_death1(4)*x(78);
        dxdt(79)=k_bind1*x(3)*(x(5)+x(8))+k_bind1*x(35)*x(2)-...
            k_bind1*x(79)*x(4)-k_intern*x(79)-k_death1(4)*x(79); % receptor-bound capBV (ThreeBac) [#/mL]
        dxdt(80)=0; % receptor-bound goiBV [#/mL]
    
        % infected by repcapBV and goiBV (TwoBac) or by repBV and goiBV (ThreeBac) (#/mL)
        dxdt(100)=k_bind1*x(2)*(x(7)+x(9))+k_bind1*x(12)*x(4)-... % receptor-bound repcapBV (TwoBac) or repBV (ThreeBac) [#/mL]
            k_bind1*x(100)*x(3)-k_intern*x(100)-k_death1(5)*x(100);
        dxdt(101)=0; % receptor-bound capBV (ThreeBac) [#/mL]
        dxdt(102)=k_bind1*x(4)*(x(5)+x(9))+k_bind1*x(58)*x(2)-...
            k_bind1*x(102)*x(3)-k_intern*x(102)-k_death1(5)*x(102); % receptor-bound goiBV [#/mL]
        
        % infected by capBV and goiBV (ThreeBac) [#/mL]
        dxdt(122)=0;  % receptor-bound repcapBV (TwoBac) or repBV (ThreeBac) [#/mL]
        dxdt(123)=k_bind1*x(3)*(x(7)+x(10))+k_bind1*x(35)*x(4)-... 
            k_bind1*x(123)*x(2)-k_intern*x(123)-k_death1(6)*x(123); % receptor-bound capBV (ThreeBac) [#/mL]
        dxdt(124)=k_bind1*x(4)*(x(6)+x(10))+k_bind1*x(58)*x(3)-...
            k_bind1*x(124)*x(2)-k_intern*x(124)-k_death1(6)*x(124); % receptor-bound goiBV [#/mL]
        
    % in cells infected by repBV, capBV, goiBV (ThreeBac) [#/mL]
    dxdt(144)=k_bind1*x(2)*(x(10)+x(11))+...
        k_bind1*x(78)*x(4)+k_bind1*x(100)*x(3)-k_intern*x(144)-...
        k_death1(7)*x(144); % receptor-bound repBV (ThreeBac) [#/mL]
    
    dxdt(145)=k_bind1*x(3)*(x(9)+x(11))+...
        k_bind1*x(79)*x(4)+k_bind1*x(123)*x(2)-k_intern*x(145)-...
        k_death1(7)*x(145); % receptor-bound capBV (ThreeBac) [#/mL]
    
    dxdt(146)=k_bind1*x(4)*(x(8)+x(11))+...
        k_bind1*x(102)*x(3)+k_bind1*x(124)*x(2)-...
        k_intern*x(146)-k_death1(7)*x(146); % receptor-bound goiBV [#/mL]
    
    %% baculovirus DNA in nucleus: transport

    % cells infected by one BV:
        % by repcapBV (TwoBac) or repBV (ThreeBac)
        dxdt(15)=-k_bind1*x(15)*(x(3)+x(4))+k_nucl*x(12)-k_death1(1)*x(15); % repcapBV DNA (TwoBac) or repBV DNA (ThreeBac) [#/mL]
        dxdt(16)=0; % capBV DNA (ThreeBac) [#/mL]
        dxdt(17)=0; % goiBV DNA  [#/mL]
    
        % by capBV (ThreeBac) 
        dxdt(37)=0;  % repcapBV DNA (TwoBac) or repBV DNA (ThreeBac) [#/mL]
        dxdt(38)=-k_bind1*x(38)*(x(2)+x(4))+k_nucl*x(35)-k_death1(2)*x(38); % capBV DNA (ThreeBac) [#/mL]
        dxdt(39)=0; % goiBV DNA  [#/mL]
    
        % by goiBV:
        dxdt(59)=0; % repcapBV DNA (TwoBac) or repBV DNA (ThreeBac) [#/mL]
        dxdt(60)=0; % capBV DNA (ThreeBac) [#/mL]
        dxdt(61)=-k_bind1*x(61)*(x(2)+x(3))+k_nucl*x(58)-k_death1(3)*x(61); % goiBV DNA  [#/mL]

    % infected by two BV 
        % infected by repBV and capBV (ThreeBac) [#/mL]
        dxdt(81)=k_bind1*x(15)*x(3)-... % repcapBV DNA (TwoBac) or repBV DNA (ThreeBac) [#/mL]
            k_bind1*x(81)*x(4)+k_nucl*x(78)-k_death1(4)*x(81);
        dxdt(82)=k_bind1*x(38)*x(2)-... % capBV DNA (ThreeBac) [#/mL]
            k_bind1*x(82)*x(4)+k_nucl*x(79)-k_death1(4)*x(82);
        dxdt(83)=0; % goiBV DNA  [#/mL]
    
        % infected by repcapBV and goiBV (TwoBac) or by repBV and goiBV (ThreeBac) (#/mL)
        dxdt(103)=k_bind1*x(15)*x(4)-... % repcapBV DNA (TwoBac) or repBV DNA (ThreeBac) [#/mL]
            k_bind1*x(103)*x(3)+k_nucl*x(100)-k_death1(5)*x(103);
        dxdt(104)=0; % capBV DNA (ThreeBac) [#/mL]
        dxdt(105)=k_bind1*x(61)*x(2)-...
            k_bind1*x(105)*x(3)+k_nucl*x(102)-k_death1(5)*x(105); % goiBV DNA  [#/mL]
    
        % infected by capBV and goiBV (ThreeBac) [#/mL]
        dxdt(125)=0; % repcapBV DNA (TwoBac) or repBV DNA (ThreeBac) [#/mL]
        dxdt(126)=k_bind1*x(38)*x(4)-... % capBV DNA (ThreeBac) [#/mL]
            k_bind1*x(126)*x(2)+k_nucl*x(123)-k_death1(6)*x(126);
        dxdt(127)=k_bind1*x(61)*x(3)-... % goiBV DNA  [#/mL]
            k_bind1*x(127)*x(2)+k_nucl*x(124)-k_death1(6)*x(127); 

    % in cells infected by repBV, capBV, goiBV (ThreeBac) [#/mL]
        dxdt(147)=k_bind1*x(81)*x(4)+k_bind1*x(103)*x(3)+k_nucl*x(144)-k_death1(7)*x(147); % repBV DNA (ThreeBac) [#/mL]
        dxdt(148)=k_bind1*x(82)*x(4)+k_bind1*x(126)*x(2)+k_nucl*x(145)-k_death1(7)*x(148); % capBV DNA (ThreeBac) [#/mL]
        dxdt(149)=k_bind1*x(105)*x(3)+k_bind1*x(127)*x(2)+k_nucl*x(146)-k_death1(7)*x(149); % goiBV DNA [#/mL]

    %% Intracellular species
    % replication, transcription and translation activation functions
    f_repl=(0.5+0.5*tanh((t-tau_repl_on)/0.3)).*(0.5+0.5*tanh((tau_repl_end-t)/1));
    f_DIE1=(0.5+0.5*tanh((t-tau_DIE1_on)/0.3)).*(0.5+0.5*tanh((tau_RNA_end-t)/1.1)); 
    f_polh=(0.5+0.5*tanh((t-tau_polh_on)/0.3)).*(0.5+0.5*tanh((tau_RNA_end-t)/1.1)); 
    f_p10=(0.5+0.5*tanh((t-tau_p10_on)/0.3)).*(0.5+0.5*tanh((tau_RNA_end-t)/1.1)); 

    for i =1:7 % i denotes the type of infected cell
        if x(i+4)>1e-16
            % viral DNA replication
            dxdt(15+(i-1)*22)=dxdt(15+(i-1)*22)+k_repl*f_repl*x(15+(i-1)*22); % nuclear repcapBV DNA (TwoBac) or repBV DNA (ThreeBac) [#/mL]
            dxdt(16+(i-1)*22)=dxdt(16+(i-1)*22)+k_repl*f_repl*x(16+(i-1)*22);  % nuclear capBV DNA (ThreeBac) [#/mL]
            dxdt(17+(i-1)*22)=dxdt(17+(i-1)*22)+k_repl*f_repl*x(17+(i-1)*22); % nuclear goiBV DNA [#/mL]
            
            % transcription 
            if BacN==3 % ThreeBac
                dxdt(18+(i-1)*22)=k_DIE1/bp_rep78*f_DIE1*x(15+(i-1)*22)-...
                    kd_rna*(x(18+(i-1)*22))-k_death1(i)*x(18+(i-1)*22); % rep78 mRNA
                dxdt(19+(i-1)*22)=k_polh/bp_rep52*f_polh*x(15+(i-1)*22)-...
                    kd_rna*(x(19+(i-1)*22))-k_death1(i)*x(19+(i-1)*22); % rep52 mRNA
                dxdt(20+(i-1)*22)=k_polh/bp_cap*f_polh*x(16+(i-1)*22)-...
                    kd_rna*x(20+(i-1)*22)-k_death1(i)*x(20+(i-1)*22); % cap mRNA
                dxdt(21+(i-1)*22)=k_p10/bp_goi_prod*f_p10*(x(17+(i-1)*22)+x(22+(i-1)*22))-...
                    kd_rna*(x(21+(i-1)*22))-k_death1(i)*x(21+(i-1)*22); % transgene mRNA
                
            else % TwoBac    
                dxdt(19+(i-1)*22)=k_polh/bp_rep78*f_polh*x(15+(i-1)*22)-...
                    kd_rna*(x(19+(i-1)*22))-k_death1(i)*x(19+(i-1)*22); % rep mRNA
                dxdt(20+(i-1)*22)=k_p10/bp_cap*f_p10*x(15+(i-1)*22)-...
                    kd_rna*x(20+(i-1)*22)-k_death1(i)*x(20+(i-1)*22); % cap mRNA
                dxdt(21+(i-1)*22)=k_p10/bp_goi_prod*f_p10*(x(17+(i-1)*22)+x(22+(i-1)*22))-...
                    kd_rna*(x(21+(i-1)*22))-k_death1(i)*x(21+(i-1)*22); % transgene mRNA
            end

            % translation 
            r(1)=k_transl/scaling*x(i+4)*(x(18+(i-1)*22))/(K_transl/scaling*x(i+4)+x(18+(i-1)*22));
            r(2)=k_transl/scaling*x(i+4)*(x(19+(i-1)*22))/(K_transl/scaling*x(i+4)+x(19+(i-1)*22));
            r(3)=k_transl/scaling*x(i+4)*(x(20+(i-1)*22))/(K_transl/scaling*x(i+4)+x(20+(i-1)*22));
            r(4)=k_transl/scaling*x(i+4)*(x(21+(i-1)*22))/(K_transl/scaling*x(i+4)+x(21+(i-1)*22));

            if BacN==2 % leaky scanning of Rep transcript in TwoBac
                r(2)=r(2)/2; 
                r(1)=r(2);
            end
            
            if sum(r)>k_transl*x(i+4)/scaling % translation machinery saturation
                r=r./sum(r)*k_transl*x(i+4)/scaling;
            end

            dxdt(23+(i-1)*22)=r(1)/bp78_transl-kd_rep*(x(23+(i-1)*22))-k_death1(i)*x(23+(i-1)*22); % rep78
            dxdt(24+(i-1)*22)=r(2)/bp52_transl-kd_rep*(x(24+(i-1)*22))-k_death1(i)*x(24+(i-1)*22); % rep52
            dxdt(25+(i-1)*22)=r(3)/bp_cap/60-k_death1(i)*x(25+(i-1)*22); % capsids
            dxdt(26+(i-1)*22)=r(4)/bp_goi_prod-kd_gfp*(x(26+(i-1)*22))-k_death1(i)*x(26+(i-1)*22); % GFP

            % transgene amplification
            dxdt(22+(i-1)*22)=k_goi_repl*2700/bp_goi/scaling*x(i+4)*...
                x(23+(i-1)*22)/(Kgoi_rep/scaling*x(i+4)+x(23+(i-1)*22))-...
                kd_dna*x(22+(i-1)*22)-k_death1(i)*x(22+(i-1)*22);

            % encapsidation
            Kencaps=min(x(22+(i-1)*22),x(25+(i-1)*22))*K_encaps_coeff+1;
            r_pack=k_pack/bp_goi*min(x(22+(i-1)*22),x(25+(i-1)*22))*...
                x(24+(i-1)*22)/(x(24+(i-1)*22)+Kencaps);

            dxdt(25+(i-1)*22)=dxdt(25+(i-1)*22)-r_pack;
            dxdt(22+(i-1)*22)=dxdt(22+(i-1)*22)-r_pack;
            dxdt(27+(i-1)*22)=r_pack-k_death1(i)*x(27+(i-1)*22);

            % nonviable cells balances
            dxdt(28+(i-1)*22) = k_death1(i)*x(4+i); % nonviable cells density
            dxdt(29+(i-1)*22) = k_death1(i)*x(22+(i-1)*22)-...
                kd_dna*x(29+(i-1)*22); % number of GOI copies in nonviable cells
            dxdt(30+(i-1)*22) = k_death1(i)*x(23+(i-1)*22)-kd_rep*(x(30+(i-1)*22)); % rep78 protein conc.
            dxdt(31+(i-1)*22) = k_death1(i)*x(24+(i-1)*22)-kd_rep*(x(31+(i-1)*22)); % rep52 protein conc.
            dxdt(32+(i-1)*22) = k_death1(i)*x(25+(i-1)*22); % empty capsids conc.
            dxdt(33+(i-1)*22) = k_death1(i)*x(27+(i-1)*22); % filled capsids conc.
        end
    end
    %% second wave
    %% infected cells  balances
    dxdt(205)=k_bind02*x(1)*x(2)-k_bind2*(x(205)*(x(3)+x(4)))-k_death2(1)*x(205); 
    dxdt(206)=k_bind02*x(1)*x(3)-k_bind2*(x(206)*(x(2)+x(4)))-k_death2(2)*x(206); 
    dxdt(207)=k_bind02*x(1)*x(4)-k_bind2*(x(207)*(x(2)+x(3)))-k_death2(3)*x(207); 
   
    dxdt(208)=k_bind2*(x(205)*x(3)+x(206)*x(2)-x(208)*x(4))-k_death2(4)*x(208); 
    dxdt(209)=k_bind2*(x(205)*x(4)+x(207)*x(2)-x(209)*x(3))-k_death2(5)*x(209); 
    dxdt(210)=k_bind2*(x(206)*x(4)+x(207)*x(3)-x(210)*x(2))-k_death2(6)*x(210); 
    
    dxdt(211)=k_bind2*(x(208)*x(4)+x(209)*x(3)+x(210)*x(2))-k_death2(7)*x(211); 
    
    %% virus bound to infected cells
    % infected by one BV
    dxdt(212)=k_bind02*x(2)*x(1)+k_bind2*x(2)*x(205)-k_bind2*x(212)*(x(3)+x(4))-...
        k_intern*x(212)-k_death2(1)*x(212);
    dxdt(213)=0;
    dxdt(214)=0;

    dxdt(234)=0;
    dxdt(235)=k_bind02*x(3)*x(1)+k_bind2*x(3)*x(206)-k_bind2*x(235)*(x(2)+x(4))-...
        k_intern*x(235)-k_death2(2)*x(235);
    dxdt(236)=0;

    dxdt(256)=0;
    dxdt(257)=0;
    dxdt(258)=k_bind02*x(4)*x(1)+k_bind2*x(4)*x(207)-k_bind2*x(258)*(x(2)+x(3))-...
        k_intern*x(258)-k_death2(3)*x(258);

    % infected by two BV 
    dxdt(278)=k_bind2*x(2)*(x(206)+x(208))+k_bind2*x(212)*x(3)-... 
        k_bind2*x(278)*x(4)-k_intern*x(278)-k_death2(4)*x(278);
    dxdt(279)=k_bind2*x(3)*(x(205)+x(208))+k_bind2*x(235)*x(2)-...
        k_bind2*x(279)*x(4)-k_intern*x(279)-k_death2(4)*x(279);
    dxdt(280)=0;

    dxdt(300)=k_bind2*x(2)*(x(207)+x(209))+k_bind2*x(212)*x(4)-... 
        k_bind2*x(300)*x(3)-k_intern*x(300)-k_death2(5)*x(300);
    dxdt(301)=0;
    dxdt(302)=k_bind2*x(4)*(x(205)+x(209))+k_bind2*x(258)*x(2)-...
        k_bind2*x(302)*x(3)-k_intern*x(302)-k_death2(5)*x(302);

    dxdt(322)=0;                                           
    dxdt(323)=k_bind2*x(3)*(x(207)+x(210))+k_bind2*x(235)*x(4)-... 
        k_bind2*x(323)*x(2)-k_intern*x(323)-k_death2(6)*x(323);
    dxdt(324)=k_bind2*x(4)*(x(206)+x(210))+k_bind2*x(258)*x(3)-...
        k_bind2*x(324)*x(2)-k_intern*x(324)-k_death2(6)*x(324);
    
    % triple infected - x(11)
    dxdt(344)=k_bind2*x(2)*(x(210)+x(211))+...
        k_bind2*x(278)*x(4)+k_bind2*x(300)*x(3)-k_intern*x(344)-k_death2(7)*x(344);
    
    dxdt(345)=k_bind2*x(3)*(x(209)+x(211))+...
        k_bind2*x(279)*x(4)+k_bind2*x(323)*x(2)-k_intern*x(345)-k_death2(7)*x(345);
    
    dxdt(346)=k_bind2*x(4)*(x(208)+x(211))+...
        k_bind2*x(302)*x(3)+k_bind2*x(324)*x(2)-k_intern*x(346)-k_death2(7)*x(346);
    
    %% virus in nucleus: transport

    % infected by one BV
    dxdt(215)=-k_bind2*x(215)*(x(3)+x(4))+k_nucl*x(212)-k_death2(1)*x(215); 
    dxdt(216)=0;
    dxdt(217)=0;

    dxdt(237)=0;
    dxdt(238)=-k_bind2*x(238)*(x(2)+x(4))+k_nucl*x(235)-k_death2(2)*x(238);
    dxdt(239)=0;

    dxdt(259)=0;
    dxdt(260)=0;
    dxdt(261)=-k_bind2*x(261)*(x(2)+x(3))+k_nucl*x(258)-k_death2(3)*x(261);

    % infected by two BV 
    dxdt(281)=k_bind2*x(215)*x(3)-... % rep_cap - x(8)
        k_bind2*x(281)*x(4)+k_nucl*x(278)-k_death2(4)*x(281);
    dxdt(282)=k_bind2*x(238)*x(2)-...
        k_bind2*x(282)*x(4)+k_nucl*x(279)-k_death2(4)*x(282);
    dxdt(283)=0;

    dxdt(303)=k_bind2*x(215)*x(4)-... % rep_goi - x(9)
        k_bind2*x(303)*x(3)+k_nucl*x(300)-k_death2(5)*x(303);
    dxdt(304)=0;
    dxdt(305)=k_bind2*x(261)*x(2)-...
        k_bind2*x(305)*x(3)+k_nucl*x(302)-k_death2(5)*x(305);

    dxdt(325)=0;
    dxdt(326)=k_bind2*x(238)*x(4)-... % cap_goi - x(10)
        k_bind2*x(326)*x(2)+k_nucl*x(323)-k_death2(6)*x(326);
    dxdt(327)=k_bind2*x(261)*x(3)-...
        k_bind2*x(327)*x(2)+k_nucl*x(324)-k_death2(6)*x(327);

    % triple infected
    dxdt(347)=k_bind2*x(281)*x(4)+k_bind2*x(303)*x(3)+k_nucl*x(344)-k_death2(7)*x(347);
    dxdt(348)=k_bind2*x(282)*x(4)+k_bind2*x(326)*x(2)+k_nucl*x(345)-k_death2(7)*x(348);
    dxdt(349)=k_bind2*x(305)*x(3)+k_bind2*x(327)*x(2)+k_nucl*x(346)-k_death2(7)*x(349);


    %% intracellular species balances
    % replication, transcription and translation activation functions
    f_repl=(0.5+0.5*tanh((t-tau_rel_on-tau_repl_on)/0.3)).*(0.5+0.5*tanh((tau_repl_end+tau_rel_on-t)/1)); % smooth transition
    f_DIE1=(0.5+0.5*tanh((t-tau_DIE1_on-tau_rel_on)/0.3)).*(0.5+0.5*tanh((tau_RNA_end+tau_rel_on-t)/1.1)); % smooth transition
    f_polh=(0.5+0.5*tanh((t-tau_polh_on-tau_rel_on)/0.3)).*(0.5+0.5*tanh((tau_RNA_end+tau_rel_on-t)/1.1)); % smooth transition
    f_p10=(0.5+0.5*tanh((t-tau_p10_on-tau_rel_on)/0.3)).*(0.5+0.5*tanh((tau_RNA_end+tau_rel_on-t)/1.1)); % smooth transition

    for i =1:7
        if x(i+204)>1e-16
            % viral DNA replication
            dxdt(215+(i-1)*22)=dxdt(215+(i-1)*22)+k_repl*f_repl*x(215+(i-1)*22); % nuclear rep DNA
            dxdt(216+(i-1)*22)=dxdt(216+(i-1)*22)+k_repl*f_repl*x(216+(i-1)*22);  % nuclear cap DNA
            dxdt(217+(i-1)*22)=dxdt(217+(i-1)*22)+k_repl*f_repl*x(217+(i-1)*22); % nuclear GOI DNA
            
            % transcription 
            if BacN==3 % ThreeBac
                dxdt(218+(i-1)*22)=k_DIE1/bp_rep78*f_DIE1*x(215+(i-1)*22)-...
                    kd_rna*(x(218+(i-1)*22))-k_death2(i)*x(218+(i-1)*22); % rep78 mRNA
                dxdt(219+(i-1)*22)=k_polh/bp_rep52*f_polh*x(215+(i-1)*22)-...
                    kd_rna*(x(219+(i-1)*22))-k_death2(i)*x(219+(i-1)*22); % rep52 mRNA
                dxdt(220+(i-1)*22)=k_polh/bp_cap*f_polh*x(216+(i-1)*22)-...
                    kd_rna*x(220+(i-1)*22)-k_death2(i)*x(220+(i-1)*22); % VP mRNA
                dxdt(221+(i-1)*22)=k_p10/bp_goi_prod*f_p10*(x(217+(i-1)*22)+x(222+(i-1)*22))-...
                    kd_rna*(x(221+(i-1)*22))-k_death2(i)*x(221+(i-1)*22); % transgene mRNA
            else % TwoBac
                dxdt(219+(i-1)*22)=k_polh/bp_rep52*f_polh*x(215+(i-1)*22)-...
                    kd_rna*(x(219+(i-1)*22))-k_death2(i)*x(219+(i-1)*22); % rep52 mRNA
                dxdt(220+(i-1)*22)=k_p10/bp_cap*f_p10*x(215+(i-1)*22)-...
                    kd_rna*x(220+(i-1)*22)-k_death2(i)*x(220+(i-1)*22); % VP mRNA
                dxdt(221+(i-1)*22)=k_p10/bp_goi_prod*f_p10*(x(217+(i-1)*22)+x(222+(i-1)*22))-...
                    kd_rna*(x(221+(i-1)*22))-k_death2(i)*x(221+(i-1)*22); % transgene mRNA
            end
    
            % translation - #/mL (scaled)
            r(1)=k_transl/scaling*x(i+204)*(x(218+(i-1)*22))/...
                (K_transl/scaling*x(i+204)+x(218+(i-1)*22));
            r(2)=k_transl/scaling*x(i+204)*(x(219+(i-1)*22))/...
                (K_transl/scaling*x(i+204)+x(219+(i-1)*22));
            r(3)=k_transl/scaling*x(i+204)*(x(220+(i-1)*22))/...
                (K_transl/scaling*x(i+204)+x(220+(i-1)*22));
            r(4)=k_transl/scaling*x(i+204)*(x(221+(i-1)*22))/...
                (K_transl/scaling*x(i+204)+x(221+(i-1)*22));

            if BacN==2
                r(2)=r(2)/2;
                r(1)=r(2);
                bp_rep78=bp_rep52;
            end
            
            if sum(r)>k_transl*x(i+204)/scaling
                r=r./sum(r)*k_transl*x(i+204)/scaling;
            end

            dxdt(223+(i-1)*22)=r(1)/bp_rep78-kd_rep*(x(223+(i-1)*22))-...
                k_death2(i)*x(223+(i-1)*22); % rep78
            dxdt(224+(i-1)*22)=r(2)/bp_rep52-kd_rep*(x(224+(i-1)*22))-...
                k_death2(i)*x(224+(i-1)*22); % rep52
            dxdt(225+(i-1)*22)=r(3)/bp_cap/60-k_death2(i)*x(225+(i-1)*22); % capsids
            dxdt(226+(i-1)*22)=r(4)/bp_goi_prod-kd_gfp*(x(226+(i-1)*22))-...
                k_death2(i)*x(226+(i-1)*22); % GFP
            
            % transgene amplification 
            dxdt(222+(i-1)*22)=k_goi_repl*2700/bp_goi/scaling*x(i+204)*...
                (x(217+(i-1)*22)+x(222+(i-1)*22))/...
                (Kgoi_dna/scaling*x(i+204)+x(217+(i-1)*22)+x(222+(i-1)*22))*...
                x(223+(i-1)*22)/...
                (Kgoi_rep/scaling*x(i+204)+x(223+(i-1)*22))-kd_dna*x(222+(i-1)*22)-...
                k_death2(i)*x(222+(i-1)*22);

            % transgene encapsidation
            Kencaps=min(x(222+(i-1)*22),x(225+(i-1)*22))*5+1;
            r_pack=k_pack/bp_goi*min(x(222+(i-1)*22),x(225+(i-1)*22))*...
                x(224+(i-1)*22)/(x(224+(i-1)*22)+Kencaps/scaling*x(i+204));

            dxdt(225+(i-1)*22)=dxdt(225+(i-1)*22)-r_pack;
            dxdt(222+(i-1)*22)=dxdt(222+(i-1)*22)-r_pack;
            dxdt(227+(i-1)*22)=r_pack-k_death2(i)*x(227+(i-1)*22); 

            % nonviable cells balances
            dxdt(228+(i-1)*22) = k_death2(i)*x(4+i); % number of nonviable cells
            dxdt(229+(i-1)*22) = k_death2(i)*x(222+(i-1)*22)-...
                kd_dna*x(229+(i-1)*22); % number of GOI copies in nonviable cells
            dxdt(230+(i-1)*22) = k_death2(i)*x(223+(i-1)*22)-kd_rep*(x(230+(i-1)*22)); % rep78 protein conc.
            dxdt(231+(i-1)*22) = k_death2(i)*x(224+(i-1)*22)-kd_rep*(x(231+(i-1)*22)); % rep52 protein conc.
            dxdt(232+(i-1)*22) = k_death2(i)*x(225+(i-1)*22); % empty capsids conc.
            dxdt(233+(i-1)*22) = k_death2(i)*x(227+(i-1)*22); % filled capsids conc.
        end
    end
end
