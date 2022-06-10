function [t,x] = infection_age2_speed

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
maxAge_goi=90;
Nage_bins_goi=maxAge_goi/DTage_goi+1;

DTage_rc=DTage_goi;
maxAge_rc=maxAge_goi;
Nage_bins_rc=maxAge_rc/DTage_rc+1;

% Simulation
x0=[C0 Vgoi0 Vrc0 0 zeros(1,Nage_bins_goi) zeros(1,Nage_bins_rc)];

k_bind=1.3e-9; % virus*mL/cell/min
mu=0.028*0; % 1/h
k_deathT=8e-5*0; % 1/h
k_deathI=0.0017*0; % 1/h
k_degrBV=7e-3; %1/h

[t,x]=ode23(@inf_model,[0:5:500],x0);

% plot(t,sum(x(:,5+Nage_bins_goi:end),2))
hold on
plot(t,x(:,4))

function dxdt = inf_model(t,x)
    dxdt=zeros(length(x),1);
    totI_goi=sum(x(5:4+Nage_bins_goi));
    totI_rc=sum(x(5+Nage_bins_goi:end));
    
    dxdt(1)=mu*x(1)-k_deathT*x(1)-k_bind*x(1)*(x(2)+x(3)); % target cells
    dxdt(2)=-k_bind*x(2)*(x(1)+totI_goi+totI_rc+x(4))-k_degrBV*x(2); % BV ITR/GOI (#/mL) 
    dxdt(3)=-k_bind*x(3)*(x(1)+totI_goi+totI_rc+x(4))-k_degrBV*x(3); % BV rep/cap (#/mL) 
    dxdt(4)=k_bind*(totI_goi*x(3)+x(2)*totI_rc)-k_deathI*x(4); % co-infected (#/mL)
    
    % infected
    dxdt(5)=k_bind*x(1)*x(2)-k_bind*x(5)*x(3)-k_deathI*x(5)-...
        x(5)/DTage_goi; % infected ITR/GOI (#/mL) 
    dxdt(4+Nage_bins_goi)=x(4+Nage_bins_goi-1)/DTage_goi-...
        k_bind*x(4+Nage_bins_goi)*x(3)-k_deathI*x(4+Nage_bins_goi);
    for i = 1:Nage_bins_goi-2
        dxdt(5+i)=(x(5+i-1)-x(5+i))/DTage_goi-k_bind*x(5+i)*...
            x(3)-k_deathI*x(5+i);
        dxdt(5+Nage_bins_goi+i)=(x(5+Nage_bins_goi+i-1)-...
            x(5+Nage_bins_goi+i))/DTage_rc-k_bind*x(5+Nage_bins_goi+i)*...
            x(2)-k_deathI*x(5+Nage_bins_goi+i);
    end 
%     dxdt(6:3+Nage_bins_goi)=(x(5:2+Nage_bins_goi)-...
%         x(6:3+Nage_bins_goi))/DTage_goi-k_bind*x(6:3+Nage_bins_goi)*...
%             x(3)-k_deathI*x(6:3+Nage_bins_goi);
% 
%     dxdt(6+Nage_bins_goi:end-1)=(x(5+Nage_bins_goi:end-2)-...
%         x(6+Nage_bins_goi:end-1))/DTage_rc-k_bind*x(6+Nage_bins_goi:end-1)*...
%         x(2)-k_deathI*x(6+Nage_bins_goi:end-1);

    dxdt(5+Nage_bins_goi)=k_bind*x(1)*x(3)-k_bind*x(5+Nage_bins_goi)*x(2)-...
        x(5+Nage_bins_goi)/DTage_rc-k_deathI*x(5+Nage_bins_goi); % infected rep/cap - age 0 (#/mL)  
    dxdt(end)=x(end-1)/DTage_rc-k_bind*x(end)*x(2)-k_deathI*x(end);
end

end