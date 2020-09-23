clear all; close all;

useWindows = 1;

alpha = 1;
gamma = 5;
tau = 2;
theta = 1;

t1 = (1/alpha)*log((exp(-alpha*tau)*(-gamma + gamma*exp(alpha*tau) + alpha*theta))/(alpha*theta));
t2 = (1/alpha)*log((exp(-alpha*tau)*theta - gamma/alpha)/(theta - gamma/alpha));
T = t1 + t2 + 2*tau;

%t1 = t1 + tau + t2;
t2 = t2 + tau + t1;

filename = 'DDE_Example1_RK3.dat';
filepath = fileparts(mfilename('fullpath'));
filepath = fullfile(filepath,filename);
if useWindows
    pathparts = strsplit(filepath,filesep);
    filepath = fullfile(pathparts{1:end-2},'Windows','Examples','DDE_Example1',pathparts{end});
end
sol = importdata(filepath);

t = sol.data(:,1);
U = sol.data(:,2:end);

plot(t,U);
xlim([t(1), t(end)]);
ylim([-0.2,4.7]);
xlabel('$t$','Interpreter','latex');
ylabel('$x(t)$','Interpreter','latex');
title(['Periodic DDE example with delay \tau = ' num2str(tau) ', and parameters: \alpha = ' num2str(alpha) ', \gamma = ' num2str(gamma) ', \theta = ' num2str(theta)]);
hold on
plot([t(1) t(end)], [theta, theta],'--', 'Color', [0.7 0.7 0.7]);
for i=0:6
    plot([t1+i*T t1+i*T], [min(U) max(U)], 'k--');
    plot([t2+i*T t2+i*T], [min(U) max(U)], 'r--');
    plot([t1+i*T t1+i*T+tau], [min(U), min(U)], 'Color', [.18, .75, .33]);
    text(t1+i*T+(0.5*tau), min(U)-0.05,'\tau', 'HorizontalAlignment', 'center','Color', [.18, .75, .33]);
    plot([t2+i*T t2+i*T+tau], [max(U), max(U)], 'Color', [.18, .75, .33]);
    text(t2+i*T+(0.5*tau), max(U)-0.05,'\tau', 'HorizontalAlignment', 'center','Color', [.18, .75, .33]);
end
legend('DRK3','\theta','t_1','t_2');