function [t,x] = infection_RK3_lagr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
% x(1) = target cell (#/mL)
% x(2) = BV ITR/GOI (#/mL)
% x(3) = BV rep/cap (#/mL)
% x(startI1:endI1) = infected ITR/GOI (#/mL)
% x(startI2:endI2) = infected rep/cap (#/mL)
% x(startCo:endCo) = co-infected cells (#/mL)

% compare new results to this: it's the first version of the lagrangian
% numerical scheme

tic
C0=1e6; % #/mL
MOIgoi=5; 
MOIrc=5;  

Vgoi0=C0*MOIgoi; % #/mL
Vrc0=C0*MOIrc; % #/mL

Dtau=2;
age_max=80;
t_plot=20;

n_bins=age_max/Dtau;
age_bins=Dtau:Dtau:age_max;
startI1=4;
endI1=3+n_bins;
startI2=4+n_bins;
endI2=3+2*n_bins;
startCo=4+2*n_bins;
endCo=3+2*n_bins+n_bins^2;

x0=[C0 Vgoi0 Vrc0 zeros(1,n_bins*2+n_bins^2)];

k_bind=1.3e-9*60; % mL/cell/h Dee and Shuler

% k_bind=1.05e-8*60; % mL/cell/h Nielsen 2003
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
                i2=2:n_bins;
                i2=i2+startCo-1;
                x(i2+n_bins*(i1-1))=x(i2+n_bins*(i1-1))+x(i2-1+n_bins*(i1-2));
    
                i1=2:n_bins-1;
                i2=n_bins;
                i2=i2+startCo-1;
                x(i2+n_bins*(i1-1))=x(i2+n_bins*(i1-1))+x(i2-1+n_bins*(i1-2));
%             
%             % reallocation
%             
                for i1=(n_bins-1):-1:2
                    for i2=(n_bins-1):-1:2
                   
                        x(startCo-1+i2+n_bins*(i1-1))=x(startCo-1+i2-1+n_bins*(i1-2));
                    end
                end

            % empty first bin
                i1=1;
                i2=1:n_bins-1;
                i2=i2+startCo-1;
                x(i2+n_bins*(i1-1))=0;

                i1=2:n_bins-1;
                i2=1;
                i2=i2+startCo-1;
                x(i2+n_bins*(i1-1))=0;

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
co_plot=reshape(co',[n_bins,n_bins]);
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

    % pre-allocate balances
    dI1dt=I1;
    dI2dt=I2;
%     dCodt=Co;
    
    % cells and virions balances
    dxdt(1)=mu*T-k_deathT*T-k_bind*T*(V1+V2); % target cells
    dxdt(2)=-k_bind*V1*C-k_degrBV*V1; % BV ITR/GOI (#/mL) 
    dxdt(3)=-k_bind*V2*C-k_degrBV*V2; % BV rep/cap (#/mL) 

    dI1dt(1)=k_bind*T*V1-k_deathI*I1(1)-k_bind*I1(1)*V2;
    dI1dt(2:end)=-k_deathI*I1(2:end)-k_bind*I1(2:end)*V2;

    dI2dt(1)=k_bind*T*V2-k_deathI*I2(1)-k_bind*I2(1)*V1;
    dI2dt(2:end)=-k_deathI*I2(2:end)-k_bind*I2(2:end)*V1;
    
    % co-infected cells
        % all bins
        dCodt=-k_deathI*Co;
        % first bin
        dCodt(1)=dCodt(1)+k_bind*(I1(1)*V2+I2(1)*V1);
        % cells freshly infected by V1
        dCodt(2:n_bins)=dCodt(2:n_bins)+k_bind*(I2(2:n_bins)*V1);
        % cells freshly infected by V2
        dCodt(1+n_bins*((2:n_bins)-1))=dCodt(1+n_bins*((2:n_bins)-1))+...
            k_bind*(I1(2:n_bins)*V2);
%     if t>0.0009 && t< 0.001
%         disp('')
%     end
    dxdt(4:end)=[dI1dt dI2dt dCodt];
      
end

end