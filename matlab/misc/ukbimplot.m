% Plots for UKBIM

close all; clear all; clc

load('interflow_line')
err_low = error;
dom_low = dom;
load('interfhigh_line')
err_high = error;
dom_high = dom;
load('interfsuperhigh_line')
err_super = error;
dom_super = dom;

indI = 80:100;
figure(1);
semilogy(1-(dom_low.reld(indI)),err_low{1}(indI),'-','DisplayName','Standard quad., 35 panels')
hold on
semilogy(1-(dom_high.reld(indI)),err_high{1}(indI),'--','DisplayName','Standard quad., 70 panels')
semilogy(1-(dom_super.reld(indI)),err_super{1}(indI),':','DisplayName','Standard quad., 140 panels')
% semilogy(1-(dom_low.reld(indI)),err_low{2}(indI),'^-','DisplayName','Special quad., 35 panels')
axis([0 0.2 1e-16 1e0]); grid on
xlabel('Relative distance to interface, $r$')
ylabel('Relative error')
legend('toggle')
publication_fig
%%
dom_low = dom;
figure(2); clf
plot(dom_low.zDrops,'-')
hold on
% plot(dom_high.z(40:4:end),'.')
axis equal

for i=1:length(dom_low.panels)
   t = dom_low.taup(dom_low.panels(i)); 
   n = -1i*t/abs(t);
   p = dom_low.tau(dom_low.panels(i));
   p1 = p-0.05*n; 
   p2 = p+0.05*n;
   line(real([p1; p2]),imag([p1; p2]),'Color','k','LineWidth',2)
end
set(gca,'visible','Off')

% plot(dom_high.tau(dom_high.panels),'x')
% % 
% plot(dom_low.tau(alpha(1:1:17)),'.','MarkerSize',25)
plot([dom_low.zDrops(1:1:16)],'.','MarkerSize',25)
axis([0.8 1.5 -0.1 0.5])

publication_fig