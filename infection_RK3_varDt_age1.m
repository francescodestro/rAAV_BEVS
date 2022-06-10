function [t,x] = infection_RK3_varDt_age1

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

Dtau=.25;

k_bind=1.3e-9*60; % mL/cell/h Dee and Shuler

% k_bind=1.05e-8*60; % mL/cell/h Nielsen 2003
mu=0.028; % 1/h
k_deathT=8e-5; % 1/h
k_deathI=0.0017; % 1/h
k_degrBV=7e-3; %1/h

% initialization
T=80;
t0=0;
hmax = 0.1* (T-t0);
rtol = 1.e-3;
atol = 1.e-6;
threshold=atol/rtol;
t=t0;
x=x0;

% first step selection
k1=inf_model(t0,[x0 0],[0 0]);
r=norm(k1./max(abs(x0),threshold), inf);
h = 0.8*rtol^(1/3)/r;

% pre-allocate arrays
xx=zeros(10e3,length(x0));
tt=zeros(10e3,1);
tt(1)=t0;
xx(1,:)=x0;
i=2;

nofailed=true;

% compute initial grid
tau1=h:min(t+h,90);
mu0=ones(1,length(tau1));
x0=[C0 Vgoi0 Vrc0 zeros(1,2) log(mu0)];

I2=zeros(1,1+length(tau1));
I2(1)=k_bind*x(2);

while t<T

    % step control
    hmin = 16*eps*abs(t);
    h=max(hmin,min(h,hmax));

    if 1.1*h >= T - t
        h = T - t;
    end

    k2=inf_model(t+h/2,x+k1*h/2);
    k3=inf_model(t+3*h/4,x+3*k2*h/4);

    t_new=t+h;
    x_new=x+h*(2*k1+3*k2+4*k3)/9;

    k4=inf_model(t_new,x_new);

    % estimate error
    e = h*(-5*k1 + 6*k2 + 8*k3 - 9*k4)/72;

    err_est=norm(e./max(max(abs(x),abs(x_new)), threshold), inf );

    if err_est <= rtol
        t=t_new;
        x=x_new;
        k1=k4;
        xx(i,:)=x;
        tt(i,:)=t;
        i=i+1;

        temp = 1.25*(err_est/rtol)^(1/3);
        h=h*min(5,1/temp);    
        nofailed=true;

        if length(tt)-i<10
            tt=[tt; zeros(1e5,1)];
            xx=[xx;zeros(1e5,length(x0))];
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
plot(t,x(:,end))

function dxdt = inf_model(t,x,I2)
    dxdt=x; % pre-allocate

    I2_tot=sum(I2);

    r_death_co=k_deathI;
    r_death_Igoi=k_deathI;
    r_death_Irc=k_deathI;

    % cells and virions balances
    dxdt(1)=mu*x(1)-k_deathT*x(1)-k_bind*x(1)*(x(2)+x(3)); % target cells
    dxdt(2)=-k_bind*x(2)*(x(1)+x(5)+I2_tot+x(4))-k_degrBV*x(2); % BV ITR/GOI (#/mL) 
    dxdt(3)=-k_bind*x(3)*(x(1)+x(5)+I2_tot+x(4))-k_degrBV*x(3); % BV rep/cap (#/mL) 
    dxdt(4)=k_bind*(x(5)*x(3)+x(2)*I2_tot)-r_death_co; % co-infected (#/mL)
    dxdt(5)=k_bind*x(1)*x(2)-r_death_Igoi-k_bind*x(5)*x(3); % infected ITR/GOI (#/mL) 
    dxdt(6)
    dxdt(7:end)=r_death_Irc+k_bind*x(2); % infected rep/cap (#/mL)
    
end

end