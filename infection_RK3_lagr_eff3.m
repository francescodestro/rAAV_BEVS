function [t,x] = infection_RK3_lagr_eff3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
% x(1) = target cell (#/mL)
% x(2) = BV ITR/GOI (#/mL)
% x(3) = BV rep/cap (#/mL)
% x(startI1:endI1) = infected ITR/GOI (#/mL)
% x(startI2:endI2) = infected rep/cap (#/mL)
% x(startCo:endCo) = co-infected cells (#/mL)

% eff3: added dependency of k_bind on infection age

tic
C0=1e6; % #/mL
MOIgoi=5; 
MOIrc=5;  

Dtau=1;
age_end_inf=20;
age_viab=80;
t_plot=60;

Vgoi0=C0*MOIgoi; % #/mL
Vrc0=C0*MOIrc; % #/mL

n_bins=age_viab/Dtau;
n_bins_inf=age_end_inf/Dtau;
age_bins=Dtau:Dtau:age_viab;

startI1=4;
endI1=3+n_bins;
startI2=4+n_bins;
endI2=3+2*n_bins;
startCo=4+2*n_bins;
endCo=3+2*n_bins+n_bins^2-(n_bins-n_bins_inf)*(n_bins-n_bins_inf+1);

% binding decay fitting parameters
tau0_bind=1.80;
beta_bind=0.148893;
k_bind0=1.05776*1e-8*60;
tauE_bind=16.10;
delta_bind=0.3321;

k_bind0=1.3e-9*60; % mL/cell/h Dee and Shuler

% k_bind=1.05e-8*60; % mL/cell/h Nielsen 2003

i=0;
ind=zeros(1,n_bins^2);
ind_sc=ind;

k_bind=ind;

% indexes co-infected cells bins and binding constant of (co-)infected cells
for i1=1:n_bins


%         n1=0.5+0.5*(tanh((current_age-tau0_bind)/.01));
%         n2=0.5+0.5*(tanh((tauE_bind-current_age)/delta_bind));
%         k_bind(i1)=k_bind0;%*((1-n1)+n1.*exp(-beta_bind*(current_age-tau0_bind))).*n2;
        
    for i2=1:n_bins


%         age_matrix(i1,i2)=max(i1,i2)*Dtau;

        if (i1-i2)^2 < n_bins_inf^2
            i=i+1;
            ind(i2+n_bins*(i1-1))=i;
%             ind_matrix(i1,i2)=i;
            ind_sc(i2+n_bins*(i1-1))=i+startCo-1;
%             ind_matrix_scaled(i1,i2)=i+startCo-1;

            current_age=max(i1,i2)*Dtau;
            n1=0.5+0.5*(tanh((current_age-tau0_bind)/.01));
            n2=0.5+0.5*(tanh((tauE_bind-current_age)/delta_bind));  
            k_bind(i2+n_bins*(i1-1))=k_bind0; %((1-n1)+n1.*exp(-beta_bind*(current_age-tau0_bind))).*n2;
        
            
        end
    end
end
k_bind=k_bind(ind>0);

% correction to reallocation when ge_end_inf = age_viability
if age_end_inf == age_viab
    a=-1;
else
    a=0;
end

x0=[C0 Vgoi0 Vrc0 zeros(1,n_bins*2+n_bins^2-(n_bins-n_bins_inf)*(n_bins-n_bins_inf+1))];

mu=0.028; % 1/h
k_deathT=8e-5; % 1/h
k_deathI=0.0017; % 1/h
k_degrBV=7e-3; %1/h

% initialization
t_final=100;
t0=0;
hmax = Dtau; % CFL condition % ode23 default: 0.1*(T-t0);
rtol = 1.e-3;
atol = 1.e-6;
threshold=atol/rtol;
t=t0;
x=x0;

% first step selection
k1=inf_model(t0,x0);
r=norm(k1./max(abs(x0),threshold), inf);
h = 0.8*rtol^(1/3)/r;

% pre-allocate arrays
xx=zeros(10e3,length(x0));
tt=zeros(10e3,1);

tt(1)=t0;
xx(1,:)=x0;
i=2;
sum_h=0;
nofailed=true;

while t<t_final

    % step control
    hmin = 16*eps*abs(t);
    h=max(hmin,min(h,hmax));

    if 1.1*h >= t_final - t
        h = t_final - t;
    elseif 1.1*h+sum_h >= Dtau
        h = Dtau - sum_h;
    end

    k2=inf_model(t+h/2,x+k1*h/2);
    k3=inf_model(t+3*h/4,x+3*k2*h/4);

    t_new=t+h;
    x_new=x+h*(2*k1+3*k2+4*k3)/9;

    k4=inf_model(t_new,x_new);

    % estimate error
    e = h*(-5*k1 + 6*k2 + 8*k3 - 9*k4)/72;
    err_est=norm(e./max(max(abs(x),abs(x_new)), threshold), inf);

    if err_est <= rtol
        t=t_new;
        x=x_new;
        k1=k4;
        xx(i,:)=x;
        tt(i,:)=t;

        i=i+1;
        sum_h=sum_h+h;

        temp = 1.25*(err_est/rtol)^(1/3);
        h=h*min(5,1/temp);    
        nofailed=true;
        
        % pre-allocation
        if length(tt)-i<10
            tt=[tt; zeros(1e5,1)];
            xx=[xx;zeros(1e5,length(x0))];
        end

        % age switch
        if sum_h >= Dtau
            sum_h=0;
            
            % reallocate I1
            x(endI1)=x(endI1)+x(endI1-1);
            x(startI1+1:endI1-1)=x(startI1:endI1-2);
            x(startI1)=0;

            % reallocate I2
            x(endI2)=x(endI2)+x(endI2-1);
            x(startI2+1:endI2-1)=x(startI2:endI2-2);
            x(startI2)=0;

            % reallocate co-infected cells
            % max age: accumulate cells
                i1=n_bins;
                i2=max(2,(n_bins-age_end_inf/Dtau+1)):n_bins;  % i2=2:n_bins;

                x(ind_sc(i2+n_bins*(i1-1)))=x(ind_sc(i2+n_bins*(i1-1)))+...
                    x(ind_sc(i2-1+n_bins*(i1-2)));
    
                i1=max(2,(n_bins-age_end_inf/Dtau+1)):n_bins-1; % i1=2:n_bins-1;
                i2=n_bins;

                x(ind_sc(i2+n_bins*(i1-1)))=x(ind_sc(i2+n_bins*(i1-1)))+...
                    x(ind_sc(i2-1+n_bins*(i1-2)));
            
             % reallocation
                for i2=(n_bins-1):-1:(age_end_inf/Dtau+1)
                    for i1=i2:-1:(i2-age_end_inf/Dtau+1)                   
                        x(ind_sc(i2+n_bins*(i1-1)))=...
                            x(ind_sc(i2-1+n_bins*(i1-2)));
                    end
                end
                for i1=(n_bins-1):-1:(age_end_inf/Dtau+1)
                    for i2=(i1-1):-1:(i1-age_end_inf/Dtau+1)                   
                        x(ind_sc(i2+n_bins*(i1-1)))=...
                            x(ind_sc(i2-1+n_bins*(i1-2)));
                    end
                end
                for i1=age_end_inf/Dtau+a:-1:2
                    for i2=(age_end_inf/Dtau)+a:-1:2
                        x(ind_sc(i2+n_bins*(i1-1)))=...
                            x(ind_sc(i2-1+n_bins*(i1-2)));
                    end
                end

            % empty first bin
                i1=1;
                i2=1:age_end_inf/Dtau; % i2=1:n_bins-1;
                x(ind_sc(i2+n_bins*(i1-1)))=0;

                i1=2:age_end_inf/Dtau;
                i2=1;
                x(ind_sc(i2+n_bins*(i1-1)))=0;

            k1=inf_model(t,x);
        end

    else
      if nofailed
        nofailed = false;
        h = max(hmin, h * max(0.5, 0.8*(rtol/err_est)^(1/3)));
      
      else
        h = max(hmin, 0.5 * h);
      end
    
    end

    if h <= hmin
        warning('Step size %e too small at t = %e.\n',h,t);
        break
    end

end
x=xx(1:i-1,:);
t=tt(1:i-1);
toc
% 
% co=x(t_plot,startCo:endCo);
co=interp1(t,x(:,startCo:endCo),t_plot);
co_plot=zeros(n_bins);
for i1=1:n_bins
    for i2=1:n_bins
        index=ind(i2+n_bins*(i1-1));
        if index >0 
            co_plot(i1,i2)=co(index);
        end
    end
end
% 
% % co=reshape(co,[n_bins,n_bins]);
bar3(co_plot)
title((t_plot))
% sum(co)

% plot(t,x(:,2))
% bar3(x(:,startI1:endI1))
% plot(t,sum(x(:,6:end),2))

function dxdt = inf_model(t,x)
    dxdt=x; % pre-allocate
    
    % species
    T=x(1);
    V1=x(2);
    V2=x(3);
    I1=x(startI1:endI1);
    I2=x(startI2:endI2);
    Co=x(startCo:endCo);

    I1_tot=sum(I1);
    I2_tot=sum(I2);
    Co_tot=sum(Co);
    C=T+I1_tot+I2_tot+Co_tot;

    
    %% uninfected cells balance
    bind1=k_bind0*T*V1; % V1 infection rate counter
    bind2=k_bind0*T*V2; % V2 infection rate counter

    dxdt(1)=mu*T-k_deathT*T-(bind1+bind2); % target cells

    %% infected cells balance

    % death and co-infection
    v1_new=k_bind(1:n_bins).*I2*V1;
    v2_new=k_bind(1:n_bins).*I1*V2;

    dI1dt=-k_deathI*I1-v2_new;
    dI2dt=-k_deathI*I2-v1_new;
   
    % new infected cells
    % virus1
    dI1dt(1)=dI1dt(1)+bind1;
    dI2dt(1)=dI2dt(1)+bind2;

    % update counters
    bind1=bind1+sum(v1_new);
    bind2=bind2+sum(v2_new);

%     % infected cells becoming co-infected
%     for ii=2:n_bins_inf
%         v1_new=k_bind(ii).*I2(ii)*V1;
%         v2_new=k_bind(ii).*I1(ii)*V2;
%         dI1dt(ii)=-k_deathI*I1(ii)-v2_new;
%         dI2dt(ii)=-k_deathI*I2(ii)-v1_new;
%         bind1=bind1+v1_new;
%         bind2=bind2+v2_new;
%     end
%         
%     % infected cells that don't bind virions anymore
%     dI1dt(n_bins_inf+1:end)=-k_deathI*I1(n_bins_inf+1:end);    
%     dI2dt(n_bins_inf+1:end)=-k_deathI*I2(n_bins_inf+1:end);
% %         disp(bind1/V1/(T+I1_tot));
%         disp(bind2/V2/(T+I2_tot));
    
    %% co-infected cells balance
    dCodt=-k_deathI*Co;
%         end
    % first bin
    dCodt(1)=dCodt(1)+k_bind(1)*(I1(1)*V2+I2(1)*V1);
    % cells freshly infected by V1
    dCodt(2:n_bins_inf)=dCodt(2:n_bins_inf)+k_bind(2:n_bins_inf).*...
        (I2(2:n_bins_inf)*V1);
    % cells freshly infected by V2
    for j=2:n_bins_inf
%             dCodt(ind(1+n_bins*((2:n_bins_inf)-1)))=dCodt(ind(1+n_bins*((2:n_bins_inf)-1)))+...
%                 k_bind*(I1(2:n_bins_inf)*V2);
        dCodt(ind(1+n_bins*(j-1)))=dCodt(ind(1+n_bins*(j-1)))+...
            k_bind(j).*(I1(j)*V2);

    end

    % virus entering already infected and co-infected cells
    % infected cells
    Co_uptake=sum(k_bind.*Co);
    bind1=bind1+V1*(sum(k_bind(1:n_bins).*I1)+Co_uptake);
    bind2=bind2+V2*(sum(k_bind(1:n_bins).*I2)+Co_uptake);
%     bind2=bind2+V2*k_bind(1)*Co(1);
%     disp(bind2/V2/(C))

        
    %% virions balance
    dxdt(2)=-bind1-k_degrBV*V1; % BV ITR/GOI (#/mL) 
    dxdt(3)=-bind2-k_degrBV*V2; % BV rep/cap (#/mL) 
%     dxdt(2)=-k_bind0*V1*C-k_degrBV*V1; % BV ITR/GOI (#/mL) 
%     dxdt(3)=-k_bind0*V2*C-k_degrBV*V2; % BV rep/cap (#/mL) 
    
    %% put together derivatives vector
    dxdt(4:end)=[dI1dt dI2dt dCodt];
%     dxdt=sparse(dxdt);
%      if t>0.0009 && t< 0.001
%         disp('')
%     end
end

end