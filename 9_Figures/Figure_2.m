%==========================================================================
% Figure_2.m
%
% Plot Figure 2 - model alpha and shoreline (extended).
% 
% Author: B. Davidson
% Last Updated: 3 October 2025
%==========================================================================

clear;
%close all;
clc;


t = 0:0.001:10;
gst = t;
g = 9.81;
s = 1/10;
da = 2;

gsT = 9.81*0.1*2;

RGB = orderedcolors("gem");
alpha1 = @(t) 0*(t<=0) + 1*da*(t>0 & t<=2) + 2*da*(t>2 & t<=4) + 3*da*(t>4 & t<=6) + 4*da*(t>6 & t<=8) + 5*da*(t>8 & t<=10);

figure(1)
clf


%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,2)
hold on
plot(t,gst,'k','linewidth',4)

plot(t,alpha1(t)-0.25*gsT,'-','linewidth',4,'color',RGB(1,:))

%Axis Labels
yticks([0 0.75*gsT 0.75*gsT+gsT])
yticklabels({'$0$','$U_{S0}$','$U_{S0} + \Delta \alpha$'})
xticks([0 2 4])
xticklabels({'$0$','$T$','$2T$'})
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',20)
xlabel('$t$','interpreter','latex')
ylabel('$\alpha$','interpreter','latex')


%Axis limits
ylim([0 5])
xlim([0 5])

legend({'\alpha_{SW}','\alpha'},'location','best')


subplot(3,2,1)
hold on

U_S0 = 0.75*gsT;
xs = U_S0*t - 0.5*g*s*t.^2;
[~, imax] = max(xs);

plot(t(1:imax),xs(1:imax),'-','linewidth',4,'color',RGB(1,:))
plot(t(imax:end),xs(imax:end),':','linewidth',4,'color',RGB(1,:))
plot(t(1:imax)+2,xs(1:imax),'-','linewidth',4,'color',RGB(1,:))
plot(t(imax:end)+2,xs(imax:end),':','linewidth',4,'color',RGB(1,:))

ylim([0 3.5])
xlim([0 4.5])

%axis labels
yticks(0:1:4)
yticklabels({})
ylabel('$x_S$','interpreter','latex')
xlabel('$t$','interpreter','latex')
xticks([0 2 4])
xticklabels({'$0$','$T$','$2T$'})
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',20)

%%%%%%%%%%%%%%%%%%%
subplot(3,2,4)
hold on
plot(t,gst,'k','linewidth',4)

plot(t,alpha1(t),'-','linewidth',4,'color',RGB(2,:))

%Axis Labels
yticks([0 gsT gsT+gsT])
yticklabels({'$0$','$U_{S0}$','$U_{S0} + \Delta \alpha$'})
xticks([0 2 4])
xticklabels({'$0$','$T$','$2T$'})
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',20)
xlabel('$t$','interpreter','latex')
ylabel('$\alpha$','interpreter','latex')

%Axis limits
ylim([0 5])
xlim([0 5])

legend({'\alpha_{SW}','\alpha'},'location','best')

subplot(3,2,3)
hold on

U_S0 = 1*gsT;
xs = U_S0*t - 0.5*g*s*t.^2;
[~, imax] = max(xs);

plot(t(1:imax),xs(1:imax),'-','linewidth',4,'color',RGB(2,:))
plot(t(imax:end),xs(imax:end),':','linewidth',4,'color',RGB(2,:))
plot(t(1:imax)+2,xs(1:imax),'-','linewidth',4,'color',RGB(2,:))
plot(t(imax:end)+2,xs(imax:end),':','linewidth',4,'color',RGB(2,:))

ylim([0 3.5])
xlim([0 4.5])

%axis labels
yticks(0:1:4)
yticklabels({})
ylabel('$x_S$','interpreter','latex')
xlabel('$t$','interpreter','latex')
xticks([0 2 4])
xticklabels({'$0$','$T$','$2T$'})
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',20)

%%%%%%%%%%%%%%%%%
subplot(3,2,6)
hold on
plot(t,gst,'k','linewidth',4)

RGB = orderedcolors("gem");
alpha1 = @(t) 0*(t<=0) + 1*da*(t>0 & t<=2) + 2*da*(t>2 & t<=4) + 3*da*(t>4 & t<=6) + 4*da*(t>6 & t<=8) + 5*da*(t>8 & t<=10);
plot(t,alpha1(t)+0.25*gsT,'-','linewidth',4,'color',RGB(3,:))

%Axis Labels
yticks([0 1.25*gsT 1.25*gsT+gsT])
yticklabels({'$0$','$U_{S0}$','$U_{S0} + \Delta \alpha$'})
xticks([0 2 4])
xticklabels({'$0$','$T$','$2T$'})
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',20)
xlabel('$t$','interpreter','latex')
ylabel('$\alpha$','interpreter','latex')

%Axis limits
ylim([0 5])
xlim([0 5])

legend({'\alpha_{SW}','\alpha'},'location','best')

subplot(3,2,5)
hold on

U_S0 = 1.25*gsT;
xs = U_S0*t - 0.5*g*s*t.^2;
[~, imax] = max(xs);

plot(t(1:imax),xs(1:imax),'-','linewidth',4,'color',RGB(3,:))
plot(t(imax:end),xs(imax:end),':','linewidth',4,'color',RGB(3,:))
plot(t(1:imax)+2,xs(1:imax),'-','linewidth',4,'color',RGB(3,:))
plot(t(imax:end)+2,xs(imax:end),':','linewidth',4,'color',RGB(3,:))

ylim([0 3.5])
xlim([0 4.5])

%axis labels
yticks(0:1:4)
yticklabels({})
ylabel('$x_S$','interpreter','latex')
xlabel('$t$','interpreter','latex')
xticks([0 2 4])
xticklabels({'$0$','$T$','$2T$'})
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',20)

% Force Painters renderer (vector-friendly)
set(gcf, 'Renderer', 'Painters');
shading flat;