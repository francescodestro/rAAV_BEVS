function [t,x] = BEVS_simulation(BacN,MOI,C0,Dt)

if BacN<2 || BacN>3
    error('The only acceptable values for BacN are 2 and 3')
end

MOIgoi=MOI(1);
MOIrep=MOI(2);
if BacN==3
    MOIcap=MOI(3);
elseif BacN==2
    MOIcap=0;
end

Vrep0=C0*MOIrep; % #/mL
Vcap0=C0*MOIcap; % #/mL
Vgoi0=C0*MOIgoi; % #/mL

x0=[C0 Vrep0 Vcap0 Vgoi0 zeros(1,361)];

options=odeset('RelTol',1e-4);

par=BEVS_parameters();

[t,x]=ode45(@BEVS_model,0:.1:Dt,x0,options,BacN,par);

scaling=par(35);
indexes=[];
for i=1:7
    indexes=[indexes [15:27 29:33]+(i-1)*22];
end
x(:,indexes)=x(:,indexes)*scaling;

end