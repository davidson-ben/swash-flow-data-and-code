clear;
close all;
clc;




%%
%load max shoreline
%compute Us
%get alpha for model

load max_shoreline.mat
g = 9.81;
s = 1/10;
h0 = 0.65;

Us = g*s*ts_max; %alpha from max shoreline time

alpha_xs = Us;

t_max_runup = sqrt((2*xs_max)/(g*s));
Us0_t_x = g*s*t_max_runup;


load alpha_measure.mat

gsT = g*s*2;
rgh = sqrt(g*h0);

alpha_gsT = alpha_measure./gsT;
alpha_rgh = alpha_measure./rgh;



%%


H = [0.1 0.125 0.15];
h0 = 0.65;


figure(1)
clf
hold on
%plot(H/h0,alpha_xs,'kx','linewidth',2,'markersize',10)
yyaxis left
plot(H/h0,alpha_gsT,'o','linewidth',2,'markersize',10)
%plot(H/h0,Us0_t_x,'k*','linewidth',2,'markersize',10)
xlabel('$H/h$ [-]','interpreter','latex')
ylabel('$\alpha/ \Delta \alpha $ [-]','interpreter','latex')
set(gca,'FontSize',25)
xlim([0.12 0.26])
xticks(0.15:0.05:0.25)
xticklabels({'0.15','0.20','0.25'})
ylim([0.9 1.1])
yticks([0.9:0.05:1.1])
yticklabels({'0.9','','1.0','','1.1',''})
set(gca,'ticklabelinterpreter','latex')

yyaxis right
plot(H/h0,alpha_rgh,'x','linewidth',2,'markersize',10)
ylim([0.7 0.8])
yticks([0.7:0.025:0.8])
yticklabels({'0.70','','0.75','','0.80'})
ylabel('$\alpha/ \sqrt{g h} $ [-]','interpreter','latex')

%legend({'\alpha_{max shoreline time}','\alpha_{measure}','\alpha_{max shoreline location}'},'location','best')

exportgraphics(gcf, "alpha_wave_plot.pdf", 'ContentType', 'vector');