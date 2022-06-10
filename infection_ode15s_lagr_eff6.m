function [t,x] = infection_ode15s_lagr_eff6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
% x(1) = target cell (#/mL)
% x(2) = BV ITR/GOI (#/mL)
% x(3) = BV rep/cap (#/mL)
% x(startI1:endI1) = infected ITR/GOI (#/mL)
% x(startI2:endI2) = infected rep/cap (#/mL)
% x(startCo:endCo) = co-infected cells (#/mL)

% eff3: added dependency of k_bind on infection age
% eff4: added binded virus balance
% eff5: added bins for virus binded to co-infected cells. bins translation
%       added, works fine
% eff6: added virus in endosome

tic
C0=1e6; % #/mL
MOIgoi=5; 
MOIrc=5;  

Dtau=1;
age_end_inf=20;
age_viab=80;
t_plot=60;
t_final=100;

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

% other state vectors: internalized virus
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

% kinetic parameters
mu=0.028; % 1/h
k_deathT=8e-5; % 1/h
k_deathI=0.0017; % 1/h
k_degrBV=7e-3; %1/h
k_bind0=1.3e-9*60; % mL/cell/h Dee and Shuler
k_endos_in1=0.023*60; %1/h
k_endos_in2=k_endos_in1; %0.023*60; %1/h
k_bind0=1.05e-8*60; % mL/cell/h Nielsen 2003

i=0;
ind=zeros(1,n_bins^2);
ind_scCo=ind;

k_bind1=ind;
ind_scB1Co=ind;
ind_scB2Co=ind;

% indexes co-infected cells bins and binding constant of (co-)infected cells
for i1=1:n_bins
        
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
            k_bind1(i2+n_bins*(i1-1))=k_bind0;%*((1-n1)+n1.*exp(-beta_bind*(current_age-tau0_bind))).*n2;
        
        end
    end
end
k_bind1=k_bind1(ind>0);
k_bind2=k_bind1;

% correction to reallocation when age_end_inf = age_viability
if age_end_inf == age_viab
    a=-1;
else
    a=0;
end

x0=[C0 Vgoi0 Vrc0 zeros(1,n_bins*2+n_bins^2-(n_bins-n_bins_inf)*(n_bins-n_bins_inf+1))];
x0_bind=zeros(1,endB2Co);
x0_endos=x0_bind;

x0=[x0 x0_bind x0_endos];

% initialization
t0=0;
t=t0;

% pre-allocate arrays
xx=zeros(10e3,length(x0));
% xx_bind=zeros(10e3,length(x0_bind));
% xx_endos=xx_bind;
tt=zeros(10e3,1);

tt(1)=t0;
xx(1,:)=x0;
xx_bind(1,:)=x0_bind;
xx_endos(1,:)=x0_endos;

i=2;
x_new=x0;

while t<t_final
        
    [t_new,x_new]=ode15s(@inf_model,[t t+Dtau],x_new);
    
    t=t_new(end);
    x=x_new(end,1:endCo);
    xb=x_new(end,endCo+1:endCo+endB2Co);
    xe=x_new(end,endB2Co+1:endB2Co+endB2Co);

    xx_bind(i,:)=xb;
    xx_endos(i,:)=xe;
    xx(i,:)=x;
    i=i+1;
  
    %% age switch
    x(endI1)=x(endI1)+x(endI1-1);
    x(startI1+1:endI1-1)=x(startI1:endI1-2);
    x(startI1)=0;
    
    % reallocate virus binded to I1
    xb(endB1)=xb(endB1)+xb(endB1-1);
    xb(startB1+1:endB1-1)=xb(startB1:endB1-2);
    xb(startB1)=0;
    
    % reallocate endosomal virus in I1
    xe(endB1)=xe(endB1)+xe(endB1-1);
    xe(startB1+1:endB1-1)=xe(startB1:endB1-2);
    xe(startB1)=0;
    
    % reallocate I2
    x(endI2)=x(endI2)+x(endI2-1);
    x(startI2+1:endI2-1)=x(startI2:endI2-2);
    x(startI2)=0;
    
    % reallocate V2 binded to I2
    xb(endB2)=xb(endB2)+xb(endB2-1);
    xb(startB2+1:endB2-1)=xb(startB2:endB2-2);
    xb(startB2)=0;
    
    xe(endB2)=xe(endB2)+xe(endB2-1);
    xe(startB2+1:endB2-1)=xe(startB2:endB2-2);
    xe(startB2)=0;
    
    % reallocate co-infected cells and virus binded to them
    % max age bins
    i1=n_bins;
    i2=max(2,(n_bins-n_bins_inf+1)):n_bins;  % i2=2:n_bins;
    
    x(ind_scCo(i2+n_bins*(i1-1)))=x(ind_scCo(i2+n_bins*(i1-1)))+...
    x(ind_scCo(i2-1+n_bins*(i1-2)));
    xb(ind_scB1Co(i2+n_bins*(i1-1)))=xb(ind_scB1Co(i2+n_bins*(i1-1)))+...
    xb(ind_scB1Co(i2-1+n_bins*(i1-2)));
    xb(ind_scB2Co(i2+n_bins*(i1-1)))=xb(ind_scB2Co(i2+n_bins*(i1-1)))+...
    xb(ind_scB2Co(i2-1+n_bins*(i1-2)));
    xe(ind_scB1Co(i2+n_bins*(i1-1)))=xe(ind_scB1Co(i2+n_bins*(i1-1)))+...
    xe(ind_scB1Co(i2-1+n_bins*(i1-2)));
    xe(ind_scB2Co(i2+n_bins*(i1-1)))=xe(ind_scB2Co(i2+n_bins*(i1-1)))+...
    xe(ind_scB2Co(i2-1+n_bins*(i1-2)));
    
    i1=max(2,(n_bins-n_bins_inf+1)):n_bins-1; % i1=2:n_bins-1;
    i2=n_bins;
    
    x(ind_scCo(i2+n_bins*(i1-1)))=x(ind_scCo(i2+n_bins*(i1-1)))+...
    x(ind_scCo(i2-1+n_bins*(i1-2)));
    xb(ind_scB1Co(i2+n_bins*(i1-1)))=xb(ind_scB1Co(i2+n_bins*(i1-1)))+...
    xb(ind_scB1Co(i2-1+n_bins*(i1-2)));
    xb(ind_scB2Co(i2+n_bins*(i1-1)))=xb(ind_scB2Co(i2+n_bins*(i1-1)))+...
    xb(ind_scB2Co(i2-1+n_bins*(i1-2)));
    xe(ind_scB1Co(i2+n_bins*(i1-1)))=xe(ind_scB1Co(i2+n_bins*(i1-1)))+...
    xe(ind_scB1Co(i2-1+n_bins*(i1-2)));
    xe(ind_scB2Co(i2+n_bins*(i1-1)))=xe(ind_scB2Co(i2+n_bins*(i1-1)))+...
    xe(ind_scB2Co(i2-1+n_bins*(i1-2)));
    
    % translation of intermediate bins
    for i2=(n_bins-1):-1:(n_bins_inf+1)
        for i1=i2:-1:(i2-age_end_inf/Dtau+1)                   
            x(ind_scCo(i2+n_bins*(i1-1)))=...
                x(ind_scCo(i2-1+n_bins*(i1-2)));
            xb(ind_scB1Co(i2+n_bins*(i1-1)))=xb(ind_scB1Co(i2-1+n_bins*(i1-2)));
            xb(ind_scB2Co(i2+n_bins*(i1-1)))=xb(ind_scB2Co(i2-1+n_bins*(i1-2)));
            xe(ind_scB1Co(i2+n_bins*(i1-1)))=xe(ind_scB1Co(i2-1+n_bins*(i1-2)));
            xe(ind_scB2Co(i2+n_bins*(i1-1)))=xe(ind_scB2Co(i2-1+n_bins*(i1-2)));
        end
    end
    for i1=(n_bins-1):-1:(n_bins_inf+1)
        for i2=(i1-1):-1:(i1-age_end_inf/Dtau+1)                   
            x(ind_scCo(i2+n_bins*(i1-1)))=...
                x(ind_scCo(i2-1+n_bins*(i1-2)));
            xb(ind_scB1Co(i2+n_bins*(i1-1)))=xb(ind_scB1Co(i2-1+n_bins*(i1-2)));
            xb(ind_scB2Co(i2+n_bins*(i1-1)))=xb(ind_scB2Co(i2-1+n_bins*(i1-2)));
            xe(ind_scB1Co(i2+n_bins*(i1-1)))=xe(ind_scB1Co(i2-1+n_bins*(i1-2)));
            xe(ind_scB2Co(i2+n_bins*(i1-1)))=xe(ind_scB2Co(i2-1+n_bins*(i1-2)));
        end
    end
    for i1=n_bins_inf+a:-1:2
        for i2=(n_bins_inf)+a:-1:2
            x(ind_scCo(i2+n_bins*(i1-1)))=...
                x(ind_scCo(i2-1+n_bins*(i1-2)));
            xb(ind_scB1Co(i2+n_bins*(i1-1)))=xb(ind_scB1Co(i2-1+n_bins*(i1-2)));
            xb(ind_scB2Co(i2+n_bins*(i1-1)))=xb(ind_scB2Co(i2-1+n_bins*(i1-2)));
            xe(ind_scB1Co(i2+n_bins*(i1-1)))=xe(ind_scB1Co(i2-1+n_bins*(i1-2)));
            xe(ind_scB2Co(i2+n_bins*(i1-1)))=xe(ind_scB2Co(i2-1+n_bins*(i1-2)));
        end
    end
    
    % reset first bin
    i1=1;
    i2=1:n_bins_inf; % i2=1:n_bins-1;
    x(ind_scCo(i2+n_bins*(i1-1)))=0;
    xb(ind_scB1Co(i2+n_bins*(i1-1)))=0;
    xb(ind_scB2Co(i2+n_bins*(i1-1)))=0;
    xe(ind_scB1Co(i2+n_bins*(i1-1)))=0;
    xe(ind_scB2Co(i2+n_bins*(i1-1)))=0;
    
    i1=2:n_bins_inf;
    i2=1;
    x(ind_scCo(i2+n_bins*(i1-1)))=0;
    xb(ind_scB1Co(i2+n_bins*(i1-1)))=0;
    xb(ind_scB2Co(i2+n_bins*(i1-1)))=0;
    xe(ind_scB1Co(i2+n_bins*(i1-1)))=0;
    xe(ind_scB2Co(i2+n_bins*(i1-1)))=0;

    x_new=[x xb xe];
   
%             % debugging: check first before and after reallocation

%             for i1=1:n_bins
%                 for i2=1:n_bins
%                     index=ind_scB1Co(i2+n_bins*(i1-1));
% %                     if index >0 
%                         x_bind_before(i1,i2)=x_bind2(index);
%                         x_bind_after(i1,i2)=x_bind(index);
%                         N_before(i1,i2)=x_bind_before(i1,i2)./x2(ind_scCo(i2+n_bins*(i1-1)));
%                         N_after(i1,i2)=x_bind_after(i1,i2)./x(ind_scCo(i2+n_bins*(i1-1)));
%     
% %                     end
%                 end
%             end
           
end


x=xx(1:i-1,:);
t=tt(1:i-1);
toc
% % 
% co=x(t_plot,startCo:endCo);

%% graphs
% if t_plot<age_viab
% % else
% co=interp1(t,x(:,startCo:endCo),t_plot);
% b1_co=interp1(t,x_endos(:,startB1Co:endB1Co),t_plot);
% b2_co=interp1(t,x_endos(:,startB2Co:endB2Co),t_plot);
    

% co=x(t_plot,startCo:endCo);
% b1_co=x_bind(t_plot,startB1Co:endB1Co);
% b2_co=x_bind(t_plot,startB2Co:endB2Co);
% 
% co_plot=zeros(n_bins);
% b1_co_plot=co_plot;
% b2_co_plot=co_plot;
% 
% for i1=1:n_bins
%     for i2=1:n_bins
%         index=ind(i2+n_bins*(i1-1));
% %         if index >0 
%             co_plot(i1,i2)=co(index);
%             b1_co_plot(i1,i2)=b1_co(index);
%             b2_co_plot(i1,i2)=b2_co(index);
% %         else
% %             co_plot(i1,i2)=0;
% %             b1_co_plot(i1,i2)=0;
% %             b2_co_plot(i1,i2)=0;
% %         end
%     end
% end

% co=reshape(co,[n_bins,n_bins]);
% bar3((b1_co_plot+b2_co_plot)./co_plot)
% U1=b1_co_plot./co_plot;
% U2=b2_co_plot./co_plot;
% % bar3(U1)
% title((t_plot))

% plot(t,x(:,2))
% bar3(x(:,startI1:endI1))
% plot(t,sum(x(:,6:end),2))

%% model
function dxdt = inf_model(t,x)
    dxdt=x; % pre-allocate
    
    % species
    T=x(1);
    V1=x(2);
    V2=x(3);
    I1=x(startI1:endI1);
    I2=x(startI2:endI2);
    Co=x(startCo:endCo);
    x_bind=x(endCo+1:endCo+endB2Co);
    x_endos=x(endCo+endB2Co+1:endCo+endB2Co+endB2Co);

    %% pre-allocate all balances 
    dI1dt=zeros(1,length(I1));
    dI2dt=dI1dt;
    dxbind_dt=zeros(1,length(x_bind));
    dxendos_dt=dxbind_dt;
    dCodt=zeros(1,length(Co)); 

    %% infection of uninfected cells
    bind1=k_bind1(1)*T*V1; % V1 infection rate counter
    bind2=k_bind2(1)*T*V2; % V2 infection rate counter

    dxdt(1)=mu*T-k_deathT*T-(bind1+bind2); % target cells balance

    dI1dt(1)=dI1dt(1)+bind1; % infected cells balance
    dI2dt(1)=dI2dt(1)+bind2; % infected cells balance

    dxbind_dt(1)=dxbind_dt(1)+bind1; % binded virus added to binded virus balance of infected cells of first age bin
    dxbind_dt(startB2)=dxbind_dt(startB2)+bind2; % binded virus added to binded virus balance of infected cells of first age bin
    
    %% I1, I2: infection, death, and intracellular species balances
    for j = 1:n_bins_inf
        % I1 and I2 death
        deathI1=k_deathI*I1(j);
        deathI2=k_deathI*I2(j);
        dI1dt(j)=dI1dt(j)-deathI1;
        dI2dt(j)=dI2dt(j)-deathI2;

        % V1, V2 binded to I1, I2: consumption for cell death and endocytosis
        Ein1=k_endos_in1*x_bind(j);
        Ein2=k_endos_in2*x_bind(startB2+j-1);
        dxbind_dt(j)=dxbind_dt(j)-Ein1-k_deathI*x_bind(j);
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)-Ein2-...
            k_deathI*x_bind(startB2+j-1); 
        dxendos_dt(j)=Ein1;
        dxendos_dt(startB2+j-1)=Ein2;

        % co-infection: virus of a different type enters infected cell
        kv1=k_bind1(j)*V1;
        kv2=k_bind2(j)*V2;
        v1_new=I2(j)*kv1; % virus 1 entering I2
        v2_new=I1(j)*kv2; % virus 2 entering I1
        bind1=bind1+v1_new; % update counters of binded virus
        bind2=bind2+v2_new; % update counters of binded virus

        dI1dt(j)=dI1dt(j)-v2_new;  % remove new co-infected cells from I1 bal
        dI2dt(j)=dI2dt(j)-v1_new;  % remove new co-infected cells from I2 bal

        dCodt(ind(1+n_bins*(j-1)))=dCodt(ind(1+n_bins*(j-1)))+v2_new; % add new co-infected cells to Co bal
        dCodt(j)=dCodt(j)+v1_new; % add new co-infected cells to Co bal

        v1_toCo=kv2*x_bind(j); % V1 binded to I1 moving to co-infected cells
        v2_toCo=kv1*x_bind(startB2+j-1); % V2 binded to I2 moving to co-infected cells

        dxbind_dt(j)=dxbind_dt(j)-v1_toCo; % move V1 binded to I1 that just got co-infected from B/I1 to B/Co balance
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)-v2_toCo;  % move V2 binded to I2 that just got co-infected from B/I2 to B/Co balance  
        
        dxbind_dt(startB1Co-1+j)=dxbind_dt(startB1Co-1+j)+v1_new; % update V1 binded to new co-infected cells
        dxbind_dt(ind_scB1Co(1+n_bins*(j-1)))=dxbind_dt(ind_scB1Co(1+n_bins*(j-1)))+... % update V1 binded to new co-infected cells
            v1_toCo;

        dxbind_dt(startB2Co-1+j)=dxbind_dt(startB2Co-1+j)+v2_toCo; % update V2 binded to new co-infected cells
        dxbind_dt(ind_scB2Co(1+n_bins*(j-1)))=dxbind_dt(ind_scB2Co(1+n_bins*(j-1)))+... % update V2 binded to new co-infected cells
            v2_new;

        % virus of the same type enters I1 and I2
        I1_uptake=k_bind1(j)*I1(j)*V1; % calculate virus uptake
        I2_uptake=k_bind2(j)*I2(j)*V2;
        dxbind_dt(j)=dxbind_dt(j)+I1_uptake; % add uptake to binded virus balance
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)+I2_uptake;
        bind1=bind1+I1_uptake; % update counters
        bind2=bind2+I2_uptake;
    end
    for j = n_bins_inf+1:n_bins
        deathI1=k_deathI*I1(j);
        deathI2=k_deathI*I2(j);

        dI1dt(j)=dI1dt(j)-deathI1;
        dI2dt(j)=dI2dt(j)-deathI2;
        
        Ein1=k_endos_in1*x_bind(j);
        Ein2=k_endos_in2*x_bind(startB2+j-1);

        dxbind_dt(j)=dxbind_dt(j)-Ein1-k_deathI*x_bind(j);
        dxbind_dt(startB2+j-1)=dxbind_dt(startB2+j-1)-Ein2-...
            k_deathI*x_bind(startB2+j-1);
        dxendos_dt(j)=Ein1;
        dxendos_dt(startB2+j-1)=Ein2;
    end
    
    % co-infected cells: death, virions binding, endocitosis
    for k=1:length(Co)
        dCodt(k)=dCodt(k)-k_deathI*Co(k);
        
        % V1 binding to co-infected cells
        Co_uptakeV1=k_bind1(k)*Co(k)*V1;
        Co_uptakeV2=k_bind2(k)*Co(k)*V2;

        Ein1=k_endos_in1*x_bind(startB1Co-1+k);
        Ein2=k_endos_in2*x_bind(startB2Co-1+k);

        dxbind_dt(startB1Co-1+k)=dxbind_dt(startB1Co-1+k)+Co_uptakeV1-...
                    k_deathI*x_bind(startB1Co-1+k)-Ein1;
        dxbind_dt(startB2Co-1+k)=dxbind_dt(startB2Co-1+k)+Co_uptakeV2-...
                    k_deathI*x_bind(startB2Co-1+k)-Ein2;

        dxendos_dt(startB1Co-1+k)=Ein1;
        dxbind_dt(startB2Co-1+k)=Ein2;
        
        bind1=bind1+Co_uptakeV1;
        bind2=bind2+Co_uptakeV2;

    end
        
    %% virions balance
    dxdt(2)=-bind1-k_degrBV*V1; % BV ITR/GOI (#/mL) 
    dxdt(3)=-bind2-k_degrBV*V2; % BV rep/cap (#/mL) 
    
    %% put together derivatives vector
    dxdt(4:end)=[dI1dt dI2dt dCodt dxbind_dt dxendos_dt];
    

end

end