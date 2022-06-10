function [t,x] = infection_age_coinf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
% x(1) = target cell (#/mL)
% x(2) = BV ITR/GOI (#/mL) 
% x(3) = BV rep/cap (#/mL) 
% x(4:4+Nage_bins_goi) = infected ITR/GOI (#/mL) 
% x(5+Nage_bins_goi:end) = infected rep/cap (#/mL)
% x(6) = co-infected (#/mL)

C0=1e6; % #/mL
MOIgoi=10; 
MOIrc=10;  

Vgoi0=C0*MOIgoi; % #/mL
Vrc0=C0*MOIrc; % #/mL

% infection age discretization grid
DTage_goi=1;
maxAge_goi=50;
Nage_bins_goi=maxAge_goi/DTage_goi+1;

DTage_rc=DTage_goi;
maxAge_rc=maxAge_goi;
Nage_bins_rc=maxAge_rc/DTage_rc+1;

% Simulation
x0=[C0 Vgoi0 Vrc0 zeros(1,Nage_bins_rc*Nage_bins_goi),...
    zeros(1,Nage_bins_goi) zeros(1,Nage_bins_rc)];

k_bind=1.3e-9; % virus*mL/cell/min
mu=0.028*0; % 1/h
k_deathT=8e-5; % 1/h
k_deathI=0.0017; % 1/h
k_degrBV=7e-3; %1/h

% indexes
start_goi=4;
end_goi=3+Nage_bins_goi;
start_rc=4+Nage_bins_goi;
end_rc=3+Nage_bins_goi+Nage_bins_rc;

[t,x]=ode23(@inf_model,[0 500],x0);

% plot(t,sum(x(:,4:3+Nage_bins_goi),2))
plot(t,sum(x(:,end_rc+1:end),2))
% hold on
% plot(t,sum(x(:,4+Nage_bins_goi:3+Nage_bins_goi+Nage_bins_rc),2))
% plot(t,x(:,end))
% plot(t,x(:,1))
% plot(t,x(:,1)+sum(x(:,4:end),2))
hold on
% plot(t,x(:,1))

function dxdt = inf_model(t,x)
    dxdt=zeros(length(x),1);
    totI_goi=sum(x(start_goi:end_goi));
    totI_rc=sum(x(start_rc:end_rc));   
%     C=reshape(,Nage_bins_goi,Nage_bins_rc);
    C=x(end_rc+1:end);
    totCo=sum(C);
%     dCdt=zeros(Nage_bins_goi,Nage_bins_rc);
    
    dxdt(1)=mu*x(1)-k_deathT*x(1)-k_bind*x(1)*(x(2)+x(3)); % target cells
    dxdt(2)=-k_bind*x(2)*(x(1)+totI_goi+totI_rc+totCo)-k_degrBV*x(2); % BV ITR/GOI (#/mL) 
    dxdt(3)=-k_bind*x(3)*(x(1)+totI_goi+totI_rc+totCo)-k_degrBV*x(3); % BV rep/cap (#/mL) 
    
    % infected ITR/GOI: first and last age bin
    dxdt(start_goi)=k_bind*x(1)*x(2)-k_bind*x(start_goi)*x(3)-...
        k_deathI*x(start_goi)-...
        x(start_goi)/DTage_goi; % infected ITR/GOI (#/mL) 
    dxdt(end_goi)=x(end_goi-1)/DTage_goi-k_bind*x(end_goi)*x(3)-...
        k_deathI*x(end_goi);
     
    % infected rep/cap: first and last age bin
    dxdt(start_rc)=k_bind*x(1)*x(3)-k_bind*x(start_rc)*x(2)-...
        x(start_rc)/DTage_rc-k_deathI*x(start_rc); % infected rep/cap - age 0 (#/mL)  
    dxdt(end_rc)=x(end_rc-1)/DTage_rc-k_bind*x(end_rc)*x(2)-...
        k_deathI*x(end_rc);

    %co-infected
    dxdt(end_rc+1)=k_bind*x(start_goi)*x(3)+k_bind*x(start_rc)*x(2)-...
        k_deathI*C(1);
    dxdt(end_rc+Nage_bins_goi)=k_bind*x(end_goi)*x(3)-k_deathI*C(Nage_bins_goi);
    dxdt(end_rc+Nage_bins_goi*Nage_bins_rc)=k_bind*x(end_rc)*x(2)-...
        k_deathI*C(Nage_bins_goi*Nage_bins_rc);

    % infected cells: intermediate age bins
    for i = 1:Nage_bins_goi-2
        dxdt(start_goi+i)=(x(start_goi+i-1)-x(4+i))/DTage_goi-...
            k_bind*x(start_goi+i)*x(3)-k_deathI*x(start_goi+i);
        dxdt(start_rc+i)=(x(start_rc+i-1)-x(start_rc+i))/...
            DTage_rc-k_bind*x(start_rc+i)*x(2)-k_deathI*x(start_rc+i);
        dxdt(end_rc+i)=k_bind*x(start_goi+i)*x(3)-k_deathI*C(end_rc+i);
        dxdt(end_rc+i+Nage_bins_rc*(Nage_bins_goi-1))=...
            k_bind*x(start_rc+i)*x(2)-...
            k_deathI*C(i+Nage_bins_rc*(Nage_bins_goi-1));
    end 

%     % co-infected
%     for i = 1 : Nage_bins_goi
%         for j = 1 : Nage_bins_rc
% 
% 
%         end
%     end
%     dxdt(end_rc+1:end)=reshape(dCdt,1,Nage_bins_goi*Nage_bins_rc);
% dxdt(end_rc+1:end_rc*2)=C(:,1);
%     dxdt(end)=k_bind*(totI_goi*x(3)+x(2)*totI_rc)-k_deathI*x(end); % co-infected (#/mL)

%          
end

end