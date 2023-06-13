% Wrapper for running the baculovirus infection and rAAV production model
clc, clear, close all

BacN=2; % 2 (TwoBac) or 3 (ThreeBac)
MOI(1)=3; % MOI: goiBV
MOI(2)=3; % MOI: repcapBV (TwoBac) or repBV (ThreeBac)
MOI(3)=3; % MOI: capBV (ThreeBac)
Dt=72; % simulation duration [h]
C0=2e6; % viable cells concentration at time of infection [#/mL]
        % hp: no nutrients limitations at time of infection
        %     the model does not consider density limitations that can occur for C0 > 5e6

[t,x] = BEVS_simulation(BacN,MOI,C0,Dt);

% output
filled=zeros(length(t),1); % filled capsids concentration [#/mL]
empty=zeros(length(t),1);  % empty capsids concentration [#/mL]
for i=1:7
    filled=filled+x(:,27+(i-1)*22)+x(:,33+(i-1)*22);
    empty=empty+x(:,25+(i-1)*22)+x(:,32+(i-1)*22);
end