function [t,x,x_bind] = infection_RK3_lagr_eff5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
% x(1) = target cell (#/mL)
% x(2) = BV ITR/GOI (#/mL)
% x(3) = BV rep/cap (#/mL)
% x(startI1:endI1) = infected ITR/GOI (#/mL)
% x(startI2:endI2) = infected rep/cap (#/mL)
% x(startCo:endCo) = co-infected bins (#/mL)

% eff3: added dependency of k_bind on infection age
% eff4: added binded virus balance

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

% first state vector
startI1=4;
endI1=3+n_bins;
startI2=4+n_bins;
endI2=3+2*n_bins;
startCo=4+2*n_bins;
endCo=3+2*n_bins+n_bins^2-(n_bins-n_bins_inf)*(n_bins-n_bins_inf+1);

% second state vector
startB1=1;
endB1=n_bins;
startB2=1+n_bins;
endB2=2*n_bins;
startB1Co=2*n_bins+1;
endB1Co=2*n_bins+n_bins^2-(n_bins-n_bins_inf)*(n_bins-n_bins_inf+1);
startB2Co=endB1Co+1;
endB2Co=endB1Co+n_bins^2-(n_bins-n_bins_inf)*(n_bins-n_bins_inf+1);

% binding decay fitting parameters
tau0_bind=1.80;
beta_bind=0.148893;
% k_bind0=1.05776*1e-8*60;
tauE_bind=16.10;
delta_bind=0.3321;

k_bind0=1.3e-9*60; % mL/cell/h Dee and Shuler
k_endos_in=0;%0.023*60; %1/h
% k_bind=1.05e-8*60; % mL/cell/h Nielsen 2003

i=0;
ind=zeros(1,n_bins^2);
ind_scCo=ind;
ind_scB1Co=ind;
ind_scB2Co=ind;

k_bind=ind;

% indexes co-infected bins bins and binding constant of (co-)infected bins
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
            ind_scCo(i2+n_bins*(i1-1))=i+startCo-1;
%             ind_matrix_scaled(i1,i2)=i+startCo-1;
            ind_scB1Co(i2+n_bins*(i1-1))=i+startB1Co-1;
            ind_scB2Co(i2+n_bins*(i1-1))=i+startB2Co-1;

            current_age=max(i1,i2)*Dtau;
            n1=0.5+0.5*(tanh((current_age-tau0_bind)/.01));
            n2=0.5+0.5*(tanh((tauE_bind-current_age)/delta_bind));  
            k_bind(i2+n_bins*(i1-1))=k_bind0; %((1-n1)+n1.*exp(-beta_bind*(current_age-tau0_bind))).*n2;
        
            
        end
    end
end
k_bind=k_bind(ind>0);

% correction to reallocation when age_end_inf = age_viability
if age_end_inf == age_viab
    a=-1;
else
    a=0;
end

x0=[C0 Vgoi0 Vrc0 zeros(1,n_bins*2+n_bins^2-(n_bins-n_bins_inf)*(n_bins-n_bins_inf+1))];
x0_bind=zeros(1,endB2Co);

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
x_bind=x0_bind;

% first step selection
[k1,k1_bind]=inf_model(t0,x0,x0_bind);
r=max(norm(k1./max(abs(x0),threshold), inf),norm(k1_bind./max(abs(x0_bind),threshold), inf));
h = 0.8*rtol^(1/3)/r;

% pre-allocate arrays
xx=zeros(10e3,length(x0));
xx_bind=zeros(10e3,length(x0_bind));
tt=zeros(10e3,1);

tt(1)=t0;
xx(1,:)=x0;
xx_bind(1,:)=x0_bind;

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

    [k2, k2_bind]=inf_model(t+h/2,x+k1*h/2,x_bind+k1_bind*h/2);
    [k3, k3_bind]=inf_model(t+3*h/4,x+3*k2*h/4,x_bind+3*k2_bind*h/4);

    t_new=t+h;
    x_new=x+h*(2*k1+3*k2+4*k3)/9;
    x_bind_new=x_bind+h*(2*k1_bind+3*k2_bind+4*k3_bind)/9;

    [k4, k4_bind]=inf_model(t_new,x_new,x_bind_new);

    % estimate error
    e = h*(-5*k1 + 6*k2 + 8*k3 - 9*k4)/72;
    e_bind=h*(-5*k1_bind + 6*k2_bind + 8*k3_bind - 9*k4_bind)/72;

    err_est=max(norm(e./max(max(abs(x),abs(x_new)), threshold), inf),...
        norm(e_bind./max(max(abs(x_bind),abs(x_bind_new)), threshold), inf));

    if err_est <= rtol
        t=t_new;
        x=x_new;
        x_bind=x_bind_new;

        k1=k4;
        k1_bind=k4_bind;

        xx(i,:)=x;
        xx_bind(i,:)=x_bind;
        tt(i,:)=t;

        i=i+1;
        sum_h=sum_h+h;

        temp = 1.25*(err_est/rtol)^(1/3);
        h=h*min(5,1/temp);    
        nofailed=true;
        
        % pre-allocate further space if necessary
        if length(tt)-i<10
            tt=[tt; zeros(1e5,1)];
            xx=[xx;zeros(1e5,length(x0))];
            xx_bind=[xx_bind;zeros(1e5,length(x0_bind))];
        end
%% age switch
        if sum_h >= Dtau
            sum_h=0;
            
            % reallocate cells and virus
            %  I1
            x(endI1)=x(endI1)+x(endI1-1);
            x(startI1+1:endI1-1)=x(startI1:endI1-2);
            x(startI1)=0;

            % reallocate virus binded to I1
            x_bind(endB1)=x_bind(endB1)+x_bind(endB1-1);
            x_bind(startB1+1:endB1-1)=x_bind(startB1:endB1-2);
            x_bind(startB1)=0;

            % reallocate I2
            x(endI2)=x(endI2)+x(endI2-1);
            x(startI2+1:endI2-1)=x(startI2:endI2-2);
            x(startI2)=0;

            % reallocate V2 binded to I2
            x_bind(endB2)=x_bind(endB2)+x_bind(endB2-1);
            x_bind(startB2+1:endB2-1)=x_bind(startB2:endB2-2);
            x_bind(startB2)=0;

            
            % reallocate co-infected cells
                i1=n_bins;
                i2=max(2,(n_bins-age_end_inf/Dtau+1)):n_bins;  % i2=2:n_bins;

                x(ind_scCo(i2+n_bins*(i1-1)))=x(ind_scCo(i2+n_bins*(i1-1)))+...
                    x(ind_scCo(i2-1+n_bins*(i1-2)));
    
                i1=max(2,(n_bins-age_end_inf/Dtau+1)):n_bins-1; % i1=2:n_bins-1;
                i2=n_bins;

                x(ind_scCo(i2+n_bins*(i1-1)))=x(ind_scCo(i2+n_bins*(i1-1)))+...
                    x(ind_scCo(i2-1+n_bins*(i1-2)));
            
             % reallocation
                for i2=(n_bins-1):-1:(age_end_inf/Dtau+1)
                    for i1=i2:-1:(i2-age_end_inf/Dtau+1)                   
                        x(ind_scCo(i2+n_bins*(i1-1)))=...
                            x(ind_scCo(i2-1+n_bins*(i1-2)));
                    end
                end
                for i1=(n_bins-1):-1:(age_end_inf/Dtau+1)
                    for i2=(i1-1):-1:(i1-age_end_inf/Dtau+1)                   
                        x(ind_scCo(i2+n_bins*(i1-1)))=...
                            x(ind_scCo(i2-1+n_bins*(i1-2)));
                    end
                end
                for i1=age_end_inf/Dtau+a:-1:2
                    for i2=(age_end_inf/Dtau)+a:-1:2
                        x(ind_scCo(i2+n_bins*(i1-1)))=...
                            x(ind_scCo(i2-1+n_bins*(i1-2)));
                    end
                end

            % empty first bin
                i1=1;
                i2=1:age_end_inf/Dtau; % i2=1:n_bins-1;
                x(ind_scCo(i2+n_bins*(i1-1)))=0;

                i1=2:age_end_inf/Dtau;
                i2=1;
                x(ind_scCo(i2+n_bins*(i1-1)))=0;

            [k1,k1_bind]=inf_model(t,x,x_bind);

%         %% age switch WITH FOR CYCLES, ATTEMPTED: reallocate bins and virus
%         if sum_h >= Dtau
%             sum_h=0;
%             
%             %% infected cells and their content
%             % max age: accumulate bins and virus
%             x(endI1)=x(endI1)+x(endI1-1);
%             x_bind(endB1)=x_bind(endB1)+x_bind(endB1-1);
%             x(endI2)=x(endI2)+x(endI2-1);
%             x_bind(endB2)=x_bind(endB2)+x_bind(endB2-1);
% 
%             % intermediate bins: translate
%             for k=n_bins-2:-1:1            
%                 x(startI1+k)=x(startI1+k-1);              
%                 x_bind(startB1+k)=x_bind(startB1+k-1);              
%                 x(startI2+k)=x(startI2+k-1);            
%                 x_bind(startB2+k)=x_bind(startB2+k-1);
%             end
% 
%             % min age: reset cell content
%             x(startI1)=0;
%             x_bind(startB1)=0;
%             x(startI2)=0;
%             x_bind(startB2)=0;
% 
% 
%             %% co-infected cells and their content
%             % max age: sum cell/virus content
%             for k=max(2,(n_bins-n_bins_inf+1)):n_bins-1
%                 i1=n_bins;
%                 i2=k;
%                 x(ind_scCo(i2+n_bins*(i1-1)))=x(ind_scCo(i2+n_bins*(i1-1)))+...
%                     x(ind_scCo(i2-1+n_bins*(i1-2)));
%                 i1=k;
%                 i2=n_bins;
%                 x(ind_scCo(i2+n_bins*(i1-1)))=x(ind_scCo(i2+n_bins*(i1-1)))+...
%                     x(ind_scCo(i2-1+n_bins*(i1-2)));
%             end
%             % update content of last bin
%             i1=n_bins;
%             i2=n_bins;
%             x(ind_scCo(i2+n_bins*(i1-1)))=x(ind_scCo(i2+n_bins*(i1-1)))+...
%                     x(ind_scCo(i2-1+n_bins*(i1-2)));
% 
%             % intermediate bins: translation
%             for k=(n_bins-1):-1:(n_bins_inf+1)
% 
%                 for kk=k-1:-1:(k-n_bins_inf+1)  
%                     i2=k;
%                     i1=kk;
%                     x(ind_scCo(i2+n_bins*(i1-1)))=...
%                             x(ind_scCo(i2-1+n_bins*(i1-2)));
%                     i1=k;
%                     i2=kk;              
%                     x(ind_scCo(i2+n_bins*(i1-1)))=...
%                         x(ind_scCo(i2-1+n_bins*(i1-2)));
%                 end
%                 i2=k;
%                 x(ind_scCo(i2+n_bins*(i1-1)))=...
%                         x(ind_scCo(i2-1+n_bins*(i1-2)));
% 
%             end
%             for i1=n_bins_inf+a:-1:2
%                 for i2=n_bins_inf+a:-1:2
%                     x(ind_scCo(i2+n_bins*(i1-1)))=...
%                         x(ind_scCo(i2-1+n_bins*(i1-2)));
%                 end
%             end
% 
%             % empty first bin
%             for k=1:n_bins_inf
% 
%                 x(ind_scCo(k))=0;
%                 x(ind_scCo(1+n_bins*(k-1)))=0;
%             end
% 
%             [k1,k1_bind]=inf_model(t,x,x_bind);
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
x_bind=xx_bind(1:i-1,:);
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

% co=reshape(co,[n_bins,n_bins]);
bar3(co_plot)
title((t_plot))
% sum(co)

% plot(t,x(:,2))
% bar3(x(:,startI1:endI1))
% plot(t,sum(x(:,6:end),2))

function [dxdt,dxbind_dt] = inf_model(t,x,x_bind)
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

    %% pre-allocate all balances 
    dI1dt=zeros(1,length(I1));
    dI2dt=dI1dt;
    dxbind_dt=zeros(1,length(x_bind));
%     dCodt=zeros(1,length(Co));
    dCodt=-k_deathI*Co;

%     dxbind_dt(end)=-k_deathI*x_bind(end);

    %% infection of uninfected bins
    bind1=k_bind0*T*V1; % V1 infection rate counter
    bind2=k_bind0*T*V2; % V2 infection rate counter

    dxdt(1)=mu*T-k_deathT*T-(bind1+bind2); % target bins balance

    dI1dt(1)=dI1dt(1)+bind1; % infected bins balance
    dI2dt(1)=dI2dt(1)+bind2; % infected bins balance

    dxbind_dt(1)=dxbind_dt(1)+bind1; % binded virus added to binded virus balance of infected bins of first age bin
    dxbind_dt(startB2)=dxbind_dt(startB2)+bind2; % binded virus added to binded virus balance of infected bins of first age bin
    
    %% infection and death of infected cells
    for j = 1:n_bins_inf
        % death
        deathI1=k_deathI*I1(j);
        deathI2=k_deathI*I2(j);
        dI1dt(j)=dI1dt(j)-deathI1;
        dI2dt(j)=dI2dt(j)-deathI2;

        dxbind_dt(j)=dxbind_dt(j)-k_endos_in*x_bind(j)-k_deathI*x_bind(j);
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)-k_endos_in*x_bind(startB2+j-1)-...
            k_deathI*x_bind(startB2+j-1); 

        % co-infection: virus of a different type enters infected cell
        kv1=k_bind(j)*V1;
        kv2=k_bind(j)*V2;

        v1_new=I2(j)*kv1; % virus 1 entering I2
        v2_new=I1(j)*kv2; % virus 2 entering I1

        dI1dt(j)=dI1dt(j)-v2_new;  % remove new co-infected cells from I1 bal
        dI2dt(j)=dI2dt(j)-v1_new;  % remove new co-infected cells from I2 bal

        dCodt(ind(1+n_bins*(j-1)))=dCodt(ind(1+n_bins*(j-1)))+... % add new co-infected cells
            v2_new;
        dCodt(ind(j))=dCodt(ind(j))+v1_new; % add new co-infected cells

        bind1=bind1+v1_new; % update counters of binded virus
        bind2=bind2+v2_new; % update counters of binded virus

        v1_toCo=kv2*x_bind(j); % virus 1 binded to I1 moving to co-infected cells
        v2_toCo=kv1*x_bind(startB2+j-1); % virus 2 binded to I2 moving to co-infected cells

        dxbind_dt(j)=dxbind_dt(j)-v1_toCo; % move binded virus from infected to new co-infected cells
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)-v2_toCo;  

        % update V1 binded to new co-infected cells
        dxbind_dt(startB1Co-1+j)=dxbind_dt(startB1Co-1+j)+v1_new;
        dxbind_dt(ind_scB1Co(1+n_bins*(j-1)))=dxbind_dt(ind_scB1Co(1+n_bins*(j-1)))+... % add new co-infected cells
            v1_toCo;

        % update V2 binded to new co-infected cells
        dxbind_dt(startB2Co-1+j)=dxbind_dt(startB2Co-1+j)+v2_new;
        dxbind_dt(ind_scB2Co(1+n_bins*(j-1)))=dxbind_dt(ind_scB1Co(1+n_bins*(j-1)))+... % add new co-infected cells
            v2_toCo;

        % virus of the same type binds infected cells
        I1_uptake=k_bind(j)*I1(j)*V1;
        I2_uptake=k_bind(j)*I2(j)*V2;
        
        dxbind_dt(j)=dxbind_dt(j)+I1_uptake;
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)+I2_uptake;
    
        bind1=bind1+I1_uptake;
        bind2=bind2+I2_uptake;

        % virus binding co-infected cells + cell death and virus endocytosis
        for v=1:n_bins
            if ind(v+n_bins*(j-1))>0
                Co_V1uptake=k_bind(v+n_bins*(j-1))*x_bind(ind_scB1Co(v+n_bins*(j-1)));
                dxbind_dt(ind_scB1Co(v+n_bins*(j-1)))=...
                    dxbind_dt(ind_scB1Co(v+n_bins*(j-1)))+Co_V1uptake-...
                    (k_endos_in+k_deathI)*x_bind(ind_scB1Co(v+n_bins*(j-1)));
    
                Co_V2uptake=k_bind(v+n_bins*(j-1))*x_bind(ind_scB2Co(v+n_bins*(j-1)));
                dxbind_dt(ind_scB2Co(v+n_bins*(j-1)))=...
                    dxbind_dt(ind_scB2Co(v+n_bins*(j-1)))+Co_V2uptake-...
                    (k_endos_in+k_deathI)*x_bind(ind_scB2Co(v+n_bins*(j-1)));
                
                bind1=bind1+Co_V1uptake;
                bind2=bind2+Co_V2uptake;
    
%                 dCodt(ind(v+n_bins*(j-1)))=dCodt(ind(v+n_bins*(j-1)))-...
%                     k_deathI*Co(ind(v+n_bins*(j-1)));
    %             dCodt(ind_scCo(v+n_bins*(j-1)))=dCodt(ind_scCo(v+n_bins*(j-1)))*...
    %                 k_deathI*Co(ind(v+n_bins*(j-1)));
            end
        end
    
    end
    % infected and co-infected cells not binding anymore virus: only death and
    % endocitosis
    for j = n_bins_inf+1:n_bins
        % I1 and I2
        deathI1=k_deathI*I1(j);
        deathI2=k_deathI*I2(j);

        dI1dt(j)=dI1dt(j)-deathI1;
        dI2dt(j)=dI2dt(j)-deathI2;

        dxbind_dt(j)=dxbind_dt(j)-k_endos_in*x_bind(j)-k_deathI*x_bind(j);
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)-k_endos_in*x_bind(startB2+j-1)-...
            k_deathI*x_bind(startB2+j-1);
        
        % Co
        for v=1:n_bins
            if ind(1+n_bins*(j-1))>0
%                 dCodt(ind(1+n_bins*(j-1)))=dCodt(ind(1+n_bins*(j-1)))-...
%                     k_deathI*Co(ind(1+n_bins*(j-1)));
            end
        % V2 binded to co-infected
%         dxbind_dt(ind_scB1Co(1+n_bins*(j-1)))=dxbind_dt(ind_scB1Co(1+n_bins*(j-1)))-... % add new co-infected cells
%             (k_endos_in+k_deathI)*x_bind(ind_scB1Co(1+n_bins*(j-1)));
% 
%         % V2 binded to co-infected
%         dxbind_dt(ind_scB2Co(1+n_bins*(j-1)))=dxbind_dt(ind_scB1Co(1+n_bins*(j-1)))-... % add new co-infected cells
%             (k_endos_in+k_deathI)*x_bind(ind_scB2Co(1+n_bins*(j-1)));
        end
    end
    
    %% virions balance
    dxdt(2)=-bind1-k_degrBV*V1; % BV ITR/GOI (#/mL) 
    dxdt(3)=-bind2-k_degrBV*V2; % BV rep/cap (#/mL) 
    
    %% put together derivatives vector
    dxdt(4:end)=[dI1dt dI2dt dCodt];

end


end