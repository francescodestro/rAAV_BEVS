function [t,x] = infection_only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
% x(1) = target cell (#/mL)
% x(2) = BV ITR/GOI (#/mL) 
% x(3) = BV rep/cap (#/mL) 
% x(4) = co-infected (#/mL)
% x(5) = infected ITR/GOI (#/mL) 
% x(6) = infected rep/cap (#/mL) 

tic
C0=1e6; % #/mL
MOIgoi=5; 
MOIrc=5;  

Vgoi0=C0*MOIgoi; % #/mL
Vrc0=C0*MOIrc; % #/mL

x0=[C0 Vgoi0 Vrc0 zeros(1,3)];

k_bind=1.3e-9*60; % mL/cell/h Dee and Shuler

% k_bind=1.05e-8*60; % mL/cell/h Nielsen 2003
mu=0.028; % 1/h
k_deathT=8e-5; % 1/h
k_deathI=0.0017; % 1/h
k_degrBV=7e-3; %1/h
k_endos_in=0.023*60; %1/h
k_endos_out=0.01*60; % 1/h
k_lys=0.01*60; % 1/h
k_nucl=0.08*60; % 1/h
k_repl=1; %1/h
k_BVrel=9.8; % virus/cell/h
tau_L=6; % h
tau_VL=18; % h
tau_VLend=80; % h
tau_repl_end=tau_VL;
k_rep_bind=5.7e-3; % molecule/cell/h
k_pack=8.18e-2; % molecule/cell/h
k_assembly=10; % fast step assumption
k_degr_rep=2.45e-2; % 1/h
k_degr_vp=2.77e-1; % 1/h

% uncertain parameters
k_AAVrel=100;
k_goi_repl=.00001;
k_rep_syn=1; % 1/h - fast step assumption
k_vp_syn=k_rep_syn;  % 1/h - fast step assumption

[t,x]=ode15s(@inf_model,0:.1:80,x0);

% plot(t,sum(x(:,4),2))
% % figure
% plot(t,x(:,7)./(x(:,4)+x(:,5)),'o')
% hold on
% plot(t,x(:,7)*C0+x(:,2),'o')
% set(gca,'xlim',[0 10])
% figure
% plot(t,x(:,7)*C0./(x(:,4)+x(:,5)))
% plot(t,x(:,4)+x(:,5))
plot(t,x(:,2))
% xlabel('Time [h]')
% ylabel('')
% set(gca,'xtick',0:10:100)
% plot(t,C0./(x(:,4)+x(:,5)))

% figure,plot(t,(x(:,1)+x(:,4)+x(:,5)+x(:,6))),hold on,plot([t(1) t(end)], [C0 C0])
toc

function dxdt = inf_model(t,x)
    dxdt=zeros(length(x),1);
%     if t>4
%         bind_decay=exp(-1/3*(t-4));
%     else
%         bind_decay=1;
%     end
%     k_bind=k_bind*bind_decay;
%     tot_cells=x(1)+x(4)+x(5)+x(6);

    % cells death rate transition (I am assuming that it switches after
    % first DNA copy enters nucleus
%     r_death_co=(k_deathT*(0.5+...
%         0.5*tanh((tau_L-t)/0.3))+k_deathI*max(1,log(x(13)+x(14)+1e-16)).*(0.5+...
%         0.5*tanh((t-tau_L)/0.3)))*x(4);
%     r_death_Igoi=(k_deathT*(0.5+...
%         0.5*tanh((tau_L-t)/0.3))+k_deathT*max(1,log(x(13)+1e-16)).*(0.5+...
%         0.5*tanh((t-tau_L)/0.3)))*x(5);
%     r_death_Irc=(k_deathT*(0.5+...
%         0.5*tanh((tau_L-t)/0.3))+k_deathT*max(1,log(x(14)+1e-16)).*(0.5+...
%         0.5*tanh((t-tau_L)/0.3)))*x(6);
%  
% 
%     r_death_co=k_deathI*0;
%     r_death_Igoi=k_deathI*0;
%     r_death_Irc=k_deathI*0;
%     k_degrBV=0; %1/h

    % cells and virions balances
    dxdt(1)=mu*x(1)-k_deathT*x(1)-k_bind*x(1)*(x(2)+x(3)); % target cells
    dxdt(2)=-k_bind*x(2)*(x(1)+x(5)+x(6)+x(4))-k_degrBV*x(2); % BV ITR/GOI (#/mL) 
    dxdt(3)=-k_bind*x(3)*(x(1)+x(5)+x(6)+x(4))-k_degrBV*x(3); % BV rep/cap (#/mL) 
    dxdt(4)=k_bind*(x(5)*x(3)+x(2)*x(6))-k_deathI*x(4); % co-infected (#/mL)
    dxdt(5)=k_bind*x(1)*x(2)-k_deathI*x(5)-k_bind*x(5)*x(3); % infected ITR/GOI (#/mL) 
    dxdt(6)=k_bind*x(1)*x(3)-k_deathI*x(6)-k_bind*x(6)*x(2); % infected rep/cap (#/mL)
    
    % avg virus bounded per cell   
%     dxdt(7)=k_bind*x(2)-k_endos_in*x(7);
%     dxdt(8)=k_bind*x(3)-k_endos_in*x(8); % BV rep/cap on surface
%     
%     % avg virus in endosome
%     dxdt(9)=k_endos_in*x(7)-k_endos_out*x(9)-k_lys*x(9);
%     dxdt(10)=k_endos_in*x(8)-k_endos_out*x(10)-k_lys*x(10);
% 
%     % avg virus in cytosol
%     dxdt(11)=k_endos_out*x(9)-k_nucl*x(11);
%     dxdt(12)=k_endos_out*x(10)-k_nucl*x(12);
% 
%     % avg virus in nucleus
%     n_L=(0.5+0.5*tanh((t-tau_L)/0.3)).*(0.5+0.5*tanh((tau_repl_end-t)/1)); % smooth transition
%     n_VL=(0.5+0.5*tanh((t-tau_VL)/1)).*(0.5+0.5*tanh((tau_VLend-t)/10)); % smooth transition
%     f_BVrepl=(1-(t-tau_L)/(tau_repl_end+3-tau_L)).*n_L; % viral DNA replication activation function
%     f_BVrel=n_VL; % BV release activation function
%     % itr/goi
%     dxdt(13)=k_nucl*x(11)+k_repl*x(13)*f_BVrepl-k_BVrel*f_BVrel*x(13)/...
%         (x(13)+x(14)+1e-12);
%     % rep/cap
%     dxdt(14)=k_nucl*x(12)+k_repl*x(14)*f_BVrepl-k_BVrel*f_BVrel*x(14)/...
%         (x(13)+x(14)+1e-12);
% 
%     % rep protein 
%     dxdt(15)=k_rep_syn*x(14)*n_VL-k_rep_bind*x(15)*x(19)-k_degr_rep*x(15)+...
%         k_pack*x(16)*x(18);
% 
%     % goi 
%     dxdt(16)=k_goi_repl*x(13)*x(15)*n_VL-k_pack*x(16)*x(18);
% 
%     % cap
%     dxdt(17)=k_vp_syn*x(14)*n_VL-k_degr_vp*x(17)-60*k_assembly*x(17);
% 
%     % CapRepComplex
%     dxdt(18)=k_rep_bind*x(15)*x(19)-k_pack*x(16)*x(18);
%     
%     % AAVempty_intra
%     dxdt(19)=k_assembly*x(17)-k_AAVrel*x(19)-k_rep_bind*x(15)*x(19);
%     
%     % AAVfilled_intra
%     dxdt(20)=k_pack*x(16)*x(18)-k_AAVrel*x(20);
% 
%     % AAVempty_EC
%     dxdt(21)=k_AAVrel*x(19)*x(4);
% 
%     % AAVfilled_EC
%     dxdt(22)=k_AAVrel*x(20)*x(4);


    
end

end