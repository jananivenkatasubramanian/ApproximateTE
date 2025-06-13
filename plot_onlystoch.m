% Plot only stochastic

blue=[0 0.4470 0.7410];
lightblue="#4DBEEE";
orange="#D95319";
red=[0.6350 0.0780 0.1840];

x=[0,0.1,0.2,0.3,0.4];

% plot(x,sims(1,3:7),':','color',lightblue,'LineWidth',1);hold on;
% plot(x,sims(2,3:7),':','color',lightblue,'LineWidth',1);
% plot(x,sims(3,3:7),':','color',lightblue,'LineWidth',1);
% plot(x,sims(4,3:7),':','color',lightblue,'LineWidth',1);
% plot(x,sims(5,3:7),':','color',lightblue,'LineWidth',1);
% plot(x,sims(6,3:7),':','color',lightblue,'LineWidth',1);
% plot(x,sims(7,3:7),':','color',lightblue,'LineWidth',1);
% plot(x,sims(8,3:7),':','color',lightblue,'LineWidth',1);
% plot(x,sims(9,3:7),':','color',lightblue,'LineWidth',1);
% h1d=plot(x,sims(10,3:7),':','color',lightblue,'LineWidth',1);
h1=plot(x,sims(36,3:7),'color',blue,'LineWidth',1.5); hold on;
plot(x,sims(36,3:7),'o','color',blue,'LineWidth',2);

% plot(x,sims(12,3:7),':','color',orange,'LineWidth',1);
% plot(x,sims(13,3:7),':','color',orange,'LineWidth',1);
% plot(x,sims(14,3:7),':','color',orange,'LineWidth',1);
% plot(x,sims(15,3:7),':','color',orange,'LineWidth',1);
% plot(x,sims(16,3:7),':','color',orange,'LineWidth',1);
% plot(x,sims(17,3:7),':','color',orange,'LineWidth',1);
% plot(x,sims(18,3:7),':','color',orange,'LineWidth',1);
% plot(x,sims(19,3:7),':','color',orange,'LineWidth',1);
% plot(x,sims(20,3:7),':','color',orange,'LineWidth',1);
% h2d=plot(x,sims(21,3:7),':','color',orange,'LineWidth',1);
h2=plot(x,sims(28,3:7),'color',red,'LineWidth',1.5);
plot(x,sims(28,3:7),'o','color',red,'LineWidth',2);

xlabel('$\omega_i$','Interpreter','latex','FontSize',13);
ylabel('Relative energy of amplitudes $a_i$','Interpreter','latex','FontSize',13);
legend([h1,h2],{'Non-stochastic',...
    'Stochastic'}, ...
    'Interpreter','latex','FontSize',13,'Location','southeast');