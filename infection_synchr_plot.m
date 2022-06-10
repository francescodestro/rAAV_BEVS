% from x to C(t,tau1,tau2)
close all

[t,x]=infection_lumped;

plot(t,x(:,13))
% [t,x] = infection;

% cells plot
plot(t,x(:,1),'linewidth',1.5)
hold on
plot(t,x(:,5),'linewidth',1.5)
plot(t,x(:,6),'linewidth',1.5)
plot(t,x(:,4),'linewidth',1.5)

xlabel('Time [h]')
ylabel('Cells concentration [#/mL])')
set(gca,'xtick',0:10:100)
set(gca,'xticklabelrotation',0,'Fontsize',18,'linewidth',1.5)
legend('Uninfected','Infected rep/cap','Infected ITR/GOI','Co-infected',...
    'location','southeast')

% virus in nucleus
figure(2)
plot(t,x(:,13),'linewidth',1.5)
hold on
plot(t,x(:,end),'linewidth',1.5)
xlabel('Time [h]')
ylabel('Virus in nucleus [#/cell])')
set(gca,'xtick',0:2:100)
set(gca,'ylim',[0 10])
set(gca,'xlim',[0 20])
set(gca,'xticklabelrotation',0,'Fontsize',18,'linewidth',1.5)
legend('ITR/GOI','rep/cap',...
    'location','southeast')
legend('rigorous','simplified kinetics',...
    'location','southeast')

return
figure
subplot(2,4,1)
plot(t,x(:,15),'linewidth',1.5)
xlabel('Time [h]')
ylabel('rep [#/cell])')
set(gca,'xticklabelrotation',0,'Fontsize',16,'linewidth',1.5)

subplot(2,4,2)
plot(t,x(:,16),'linewidth',1.5)
xlabel('Time [h]')
ylabel('goi copies [#/cell])')
set(gca,'xticklabelrotation',0,'Fontsize',16,'linewidth',1.5)

subplot(2,4,3)
plot(t,x(:,17),'linewidth',1.5)
xlabel('Time [h]')
ylabel('vp [#/cell])')
set(gca,'xticklabelrotation',0,'Fontsize',16,'linewidth',1.5)

subplot(2,4,4)
plot(t,x(:,18),'linewidth',1.5)
xlabel('Time [h]')
ylabel('CapRepComplex [#/cell])')
set(gca,'xticklabelrotation',0,'Fontsize',16,'linewidth',1.5)

subplot(2,4,5)
plot(t,x(:,19),'linewidth',1.5)
xlabel('Time [h]')
ylabel('IC AAVempty [#/cell])')
set(gca,'xticklabelrotation',0,'Fontsize',16,'linewidth',1.5)

subplot(2,4,6)
plot(t,x(:,20),'linewidth',1.5)
xlabel('Time [h]')
ylabel('IC AAVfilled [#/cell])')
set(gca,'xticklabelrotation',0,'Fontsize',16,'linewidth',1.5)

subplot(2,4,7)
plot(t,x(:,21),'linewidth',1.5)
xlabel('Time [h]')
ylabel('EC AAVempty [#/ml])')
set(gca,'xticklabelrotation',0,'Fontsize',16,'linewidth',1.5)
    
subplot(2,4,8)
plot(t,x(:,22),'linewidth',1.5)
xlabel('Time [h]')
ylabel('EC AAVfilled [#/ml])')
set(gca,'xticklabelrotation',0,'Fontsize',16,'linewidth',1.5)

figure
subplot(2,2,1)
plot(t,x(:,7),'linewidth',1.5)
xlabel('Time [h]')
ylabel('Surface-binded BV [#/cell])')
set(gca,'xtick',0:10:100)
set(gca,'xticklabelrotation',0,'Fontsize',18,'linewidth',1.5)

subplot(2,2,2)
plot(t,x(:,9),'linewidth',1.5)
xlabel('Time [h]')
ylabel('Endosomal BV [#/cell])')
set(gca,'xtick',0:10:100)
set(gca,'xticklabelrotation',0,'Fontsize',18,'linewidth',1.5)

subplot(2,2,3)
plot(t,x(:,11),'linewidth',1.5)
xlabel('Time [h]')
ylabel('Cytosolic BV [#/cell])')
set(gca,'xtick',0:10:100)
set(gca,'xticklabelrotation',0,'Fontsize',18,'linewidth',1.5)

subplot(2,2,4)
plot(t,x(:,13),'linewidth',1.5)
hold on
plot(t,x(:,end),'linewidth',1.5)
xlabel('Time [h]')
ylabel('Virus in nucleus [#/cell])')
set(gca,'xtick',0:10:100)
set(gca,'xticklabelrotation',0,'Fontsize',18,'linewidth',1.5)
