function [t,x] = infection_age1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
% x(1) = target cell (#/mL)
% x(2) = BV ITR/GOI (#/mL) 
% x(3) = BV rep/cap (#/mL) 
% x(4) = co-infected (#/mL)
% x(5) = infected ITR/GOI (#/mL) 
% x(6) = infected rep/cap (#/mL) 

C0=1e6; % #/mL
MOIgoi=5; 
MOIrc=5;  

Vgoi0=C0*MOIgoi; % #/mL
Vrc0=C0*MOIrc; % #/mL

% infection age discretization grid
DTage=1;
Agemax=90;
Nage_bins=Agemax/DTage+1;

% Simulation
x0=[C0 Vgoi0 Vrc0 0 0 zeros(1,Nage_bins)];

k_bind=1.3e-9*60; % virus*mL/cell/min
mu=0.028*0; % 1/h
k_deathT=8e-5*0; % 1/h
k_deathI=1.7e-3*0; % 1/h
k_degrBV=7e-3*0; %1/h

[t,x]=ode23(@inf_model,[0 500],x0);
plot(x(:,))
plot(t,sum(x(:,6:end),2))

function dxdt = inf_model(t,x)
    dxdt=zeros(length(x),1);
    totI_rc=sum(x(6:end));
    dxdt(1)=mu*x(1)-k_deathT*x(1)-k_bind*x(1)*(x(2)+x(3)); % target cells
    dxdt(2)=-k_bind*x(2)*(x(1)+x(5)+totI_rc+x(4))-k_degrBV*x(2); % BV ITR/GOI (#/mL) 
    dxdt(3)=-k_bind*x(3)*(x(1)+x(5)+totI_rc+x(4))-k_degrBV*x(3); % BV rep/cap (#/mL) 
    dxdt(4)=k_bind*(x(5)*x(3)+x(2)*totI_rc)-k_deathI*x(4); % co-infected (#/mL)
    dxdt(5)=k_bind*x(1)*x(2)-k_bind*x(5)*x(3)-k_deathI*x(5); % infected ITR/GOI (#/mL) 
    dxdt(6)=k_bind*x(1)*x(3)-k_bind*x(6)*x(2)-x(6)/DTage-k_deathI*x(6); % infected rep/cap - age 0 (#/mL)  
    for i = 1:Nage_bins-2
        dxdt(6+i)=(x(6+i-1)-x(6+i))/DTage-k_bind*x(6+i)*x(2)-k_deathI*x(6+i);
    end
%     dxdt(7:end-1)=(x(6:end-2)-x(7:end-1))/DTage-k_bind*x(7:end-1)*x(3)-k_deathI*x(7:end-1);
    dxdt(end)=(x(end-1))/DTage-k_bind*x(end)*x(2)-k_deathI*x(end);
end

end