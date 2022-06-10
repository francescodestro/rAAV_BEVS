function [t,x,h] = infection_age_coinf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
% x(1) = target cell (#/mL)
% x(2) = BV ITR/GOI (#/mL) 
% x(3) = BV rep/cap (#/mL) 
% x(4:4+Nage_bins_goi) = infected ITR/GOI (#/mL) 
% x(5+Nage_bins_goi:end) = infected rep/cap (#/mL)
% x(end_rc+1:end) = co-infected (#/mL)

C0=1e6; % #/mL
MOIgoi=5; 
MOIrc=5;  

Vgoi0=C0*MOIgoi; % #/mL
Vrc0=C0*MOIrc; % #/mL

% infection age discretization grid
DTage_goi=1;
maxAge_goi=100;
Nage_bins_goi=maxAge_goi/DTage_goi+1;

DTage_rc=DTage_goi;
maxAge_rc=maxAge_goi;
Nage_bins_rc=maxAge_rc/DTage_rc+1;

age_intervals=0:DTage_goi:maxAge_goi;
age_labels=[];
for i=1:Nage_bins_rc-1
    age_labels=[age_labels {[num2str(age_intervals(i)) '-' ,...
        num2str(age_intervals(i+1))]}];
end
age_labels(end+1)={['>' num2str(maxAge_rc)]};

% Simulation
x0=[C0 Vgoi0 Vrc0 zeros(1,Nage_bins_rc*Nage_bins_goi),...
    zeros(1,Nage_bins_goi) zeros(1,Nage_bins_rc)]';

k_bind=1.3e-9*60; % virus*mL/cell/h
mu=0.028; % 1/h
k_deathT=8e-5; % 1/h
k_deathI=0.0017; % 1/h
k_degrBV=7e-3; %1/h

% indexes
start_goi=4;
end_goi=3+Nage_bins_goi;
start_rc=4+Nage_bins_goi;
end_rc=3+Nage_bins_goi+Nage_bins_rc;

% inf_model(0,x0)
tic
[t,x]=ode23(@inf_model,0:1:100,x0);
toc
% figure
% plot(t,sum(x(:,4:3+Nage_bins_goi),2))
% figure(1)
% plot(t,sum(x(:,end_rc+1:end),2))
% xlabel('Time [h]')
% ylabel('Total co-infected [#/mL]')

% hold on
% plot(t,sum(x(:,4+Nage_bins_goi:3+Nage_bins_goi+Nage_bins_rc),2))
% plot(t,x(:,end))
% plot(t,x(:,1))
% plot(t,x(:,1)+sum(x(:,4:end),2))
% hold on
% % plot(t,x(:,1))
% 
% figure(2)
% plot(t,x(:,start_goi:end_goi))
% xlabel('Time [h]')
% ylabel('Infected ITR/GOI BV [#/mL]')
% legend(age_labels)
% 
% figure(3)
% plot(t,x(:,start_rc:end_rc))
% xlabel('Time [h]')
% ylabel('Infected rep/cap BV [#/mL]')
% legend(age_labels)
% 
% 
% figure(4)
% tt=30;
% tau=x(tt,end_rc+1:end);
% tau2=reshape(tau',Nage_bins_goi,Nage_bins_rc);
% bar3(tau2)
% title(['t = ' num2str(t(tt)) 'h'])
% set(gca,'xticklabel',age_labels,'yticklabel',age_labels)
% xlabel('Infection age, ITR/GOI BV [h]')
% ylabel('Infection age, rep/cap BV [h]')
% zlabel('Co-infected cells [#/mL]')
% 
figure
tt=50;
tau=x(tt,end_rc+1:end);
tau2=reshape(tau',Nage_bins_goi,Nage_bins_rc);

% step=1/DTage_rc;

% tau2=tau2(1:step:end,1:step:end);

h=bar3(tau2);
% for i = 1:length(h)
%     h(i).LineWidth = .01;
%     h(i).EdgeColor = 'none';
%     h(i).MarkerSize=.1;
% end
title(['t = ' num2str(t(tt)) 'h'])
% set(gca,'xticklabel',age_labels,'yticklabel',age_labels,'linewidth',1.5)
% xlabel('Infection age, ITR/GOI BV [h]')
% ylabel('Infection age, rep/cap BV [h]')
xlabel('t_{inf1} [hpi]')
ylabel('t_{inf2} [hpi]')
zlabel('Co-infected cells [#/mL]')
set(gca,'linewidth',1.5,'fontsize',18)

% figure(5)

function dxdt = inf_model(t,x)
    dxdt=zeros(length(x),1);
    totI_goi=sum(x(start_goi:end_goi));
    totI_rc=sum(x(start_rc:end_rc));   
    C=x(end_rc+1:end);
    totCo=sum(C);
    
    % target cells and virions balances
    dxdt(1)=mu*x(1)-k_deathT*x(1)-k_bind*x(1)*(x(2)+x(3)); % target cells
    dxdt(2)=-k_bind*x(2)*(x(1)+totI_goi+totI_rc+totCo)-k_degrBV*x(2); % BV ITR/GOI (#/mL) 
    dxdt(3)=-k_bind*x(3)*(x(1)+totI_goi+totI_rc+totCo)-k_degrBV*x(3); % BV rep/cap (#/mL) 
    
    % infected ITR/GOI: first age bin
    dxdt(start_goi)=k_bind*x(1)*x(2)-k_bind*x(start_goi)*x(3)-...
        k_deathI*x(start_goi)-...
        x(start_goi)/DTage_goi; % infected ITR/GOI (#/mL) 
    
    % infected rep/cap: first age bin
    dxdt(start_rc)=k_bind*x(1)*x(3)-k_bind*x(start_rc)*x(2)-...
        x(start_rc)/DTage_rc-k_deathI*x(start_rc); % infected rep/cap - age 0 (#/mL)  

    % infected cells: intermediate age bins
    for i = 1:Nage_bins_goi-1
        % infected by only one BV
        dxdt(start_goi+i)=(x(start_goi+i-1)-x(start_goi+i))/DTage_goi-...
            k_bind*x(start_goi+i)*x(3)-k_deathI*x(start_goi+i);
        dxdt(start_rc+i)=(x(start_rc+i-1)-x(start_rc+i))/...
            DTage_rc-k_bind*x(start_rc+i)*x(2)-k_deathI*x(start_rc+i);
        
%         % co-infected
%         % first row: newly co-infected by rep/cap BV
%         dxdt(end_rc+i+1)=k_bind*x(start_goi+i)*x(3)-k_deathI*C(1+i)-...
%             C(1+i)/DTage_rc;
%         disp(i+1)
%         % first column: newly co-infected by ITR/GOI BV
%         dxdt(end_rc+1+i*Nage_bins_goi)=...
%             k_bind*x(start_rc+i)*x(2)-...
%             k_deathI*C(1+i*Nage_bins_goi)-C(1+i*Nage_bins_goi)/DTage_rc;
%         disp(1+i*Nage_bins_goi)
%         % middle cells: infection age increases
%         dxdt(end_rc+(i+1)+Nage_bins_goi*(1:Nage_bins_rc-1))=...
%             (C(i+Nage_bins_goi*(0:Nage_bins_rc-2))-...
%             C(i+1+Nage_bins_goi*(1:Nage_bins_rc-1)))/DTage_goi-...
%             -k_deathI*C(i+1+Nage_bins_goi*(1:Nage_bins_rc-1));
%         disp((i+1)+Nage_bins_goi*(1:Nage_bins_rc-1))

    end 

    % infected: last age bin correction to eliminate outlet component
    dxdt(end_rc)=dxdt(end_rc)+x(end_rc)/DTage_rc;
    dxdt(end_goi)=dxdt(end_goi)+x(end_goi)/DTage_goi;

%     % co-infected: last column correction to eliminate outlet component
%     dxdt(end_rc+Nage_bins_rc+Nage_bins_rc*(0:Nage_bins_rc-1))=...
%         dxdt(end_rc+Nage_bins_rc+Nage_bins_rc*(0:Nage_bins_rc-1))+...
%         x(end_rc+Nage_bins_rc+Nage_bins_rc*(0:Nage_bins_rc-1))/DTage_goi;
%     disp(Nage_bins_rc+Nage_bins_rc*(0:Nage_bins_rc-1))
%     % co-infected: last row correction to eliminate outlet component
%     dxdt(end_rc+(1:Nage_bins_rc-1)+Nage_bins_rc*(Nage_bins_rc-1))=...
%         dxdt(end_rc+(1:Nage_bins_rc-1)+Nage_bins_rc*(Nage_bins_rc-1))+...
%         x(end_rc+(1:Nage_bins_rc-1)+Nage_bins_rc*(Nage_bins_rc-1))/DTage_goi;
%     disp((1:Nage_bins_rc-1)+Nage_bins_rc*(Nage_bins_rc-1))


    %% co-infected
%     if t>10 && t>10.001
%         disp('qui')
%     end
    % first bin
    dxdt(end_rc+1)=k_bind*x(start_goi)*x(3)+k_bind*x(start_rc)*x(2)-...
        k_deathI*C(1)-C(1)/DTage_goi;
%     dxdt(end)=dxdt(end)+C(1)/DTage_goi-k_deathI*C(end);
    % middle cells
    for i = 2:Nage_bins_goi-1
        for j = 2:Nage_bins_goi-1
            dxdt(end_rc+j+(i-1)*Nage_bins_goi)=(x(end_rc+j-1+(i-2)*Nage_bins_goi)-...
                x(end_rc+j+(i-1)*Nage_bins_goi))/DTage_goi-...
                k_deathI*x(end_rc+j+(i-1)*Nage_bins_goi);
%             disp(j+(i-1)*Nage_bins_goi)
        end
        j=Nage_bins_goi;
        dxdt(end_rc+j+(i-1)*Nage_bins_goi)=x(end_rc+j-1+(i-2)*Nage_bins_goi)/...
                DTage_goi-k_deathI*x(end_rc+j+(i-1)*Nage_bins_goi);
%         disp(j+(i-1)*Nage_bins_goi)
    end
    i=Nage_bins_goi;
    for j=2:Nage_bins_goi
        dxdt(end_rc+j+(i-1)*Nage_bins_goi)=x(end_rc+j-1+(i-2)*Nage_bins_goi)/...
                 DTage_goi-k_deathI*x(end_rc+j+(i-1)*Nage_bins_goi);
%         disp(j+(i-1)*Nage_bins_goi)
    end

    % first row
    i=1;
    for j=2:Nage_bins_goi-1
        dxdt(end_rc+j+(i-1)*Nage_bins_goi)=k_bind*x(3)*x(start_goi+j-1)-...
            k_deathI*x(end_rc+j+(i-1)*Nage_bins_goi)-x(end_rc+j+(i-1)*Nage_bins_goi)/DTage_goi;
%         disp(j+(i-1)*Nage_bins_goi)
    end

    j=Nage_bins_goi;
    dxdt(end_rc+j+(i-1)*Nage_bins_goi)=k_bind*x(3)*x(start_goi+j-1)-...
            k_deathI*x(end_rc+j+(i-1)*Nage_bins_goi);
    
    % first column
    j=1;
    for i=2:Nage_bins_goi-1
        dxdt(end_rc+j+(i-1)*Nage_bins_goi)=k_bind*x(2)*x(start_rc+i-1)-...
            k_deathI*x(end_rc+j+(i-1)*Nage_bins_goi)-x(end_rc+j+(i-1)*Nage_bins_goi)/DTage_goi;
%         disp(j+(i-1)*Nage_bins_goi)
    end
    i=Nage_bins_goi;
    dxdt(end_rc+j+(i-1)*Nage_bins_goi)=k_bind*x(2)*x(start_rc+i-1)-...
            k_deathI*x(end_rc+j+(i-1)*Nage_bins_goi);



end

end