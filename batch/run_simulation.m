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

% output - look at readme.md for detailed information on output

filled=zeros(length(t),1); % filled capsids concentration (in viable+nonviable cells) [#/mL]
empty=zeros(length(t),1);  % empty capsids concentration (in viable+nonviable cells) [#/mL]
rep78=zeros(length(t),1);  % rep78 concentration (in viable+nonviable cells) [#/mL]
rep52=zeros(length(t),1);  % rep78 concentration (in viable+nonviable cells) [#/mL]
vg=zeros(length(t),1);     % non-encapsidated vector genome concentration (in viable+nonviable cells) [#/mL]

T=x(:,1); % uninfected cells concentration

for i=1:7
    filled=filled+x(:,27+(i-1)*22)+x(:,33+(i-1)*22);
    empty=empty+x(:,25+(i-1)*22)+x(:,32+(i-1)*22);
    rep78=rep78+x(:,23+(i-1)*22)+x(:,30+(i-1)*22);
    rep52=rep52+x(:,24+(i-1)*22)+x(:,31+(i-1)*22);
    vg=vg+x(:,22+(i-1)*22)+x(:,29+(i-1)*22);
end