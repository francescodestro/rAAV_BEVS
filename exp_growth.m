% exponential growth
close all
hold on

t=6:.1:20;
tau=5;
k=1/tau;

for y0=[.1 1 5]
    y=y0*exp(k*(t-6));
    
    plot(t,y,'linewidth',1.5)
end
set(gca,'xticklabelrotation',0,'fontsize',18,'linewidth',1.5)
xlabel('Time [hpi]')
ylabel('DNA copies [#/cell]')
box on

legend('0.1','1','5','location','northwest')

