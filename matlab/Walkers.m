clear all;
close all;
figure(1);
walkers = load('walkers.txt');
ln = size(walkers);
n = ln(1);
l = ln(2);
x = 0.1:0.1:n*0.1;

hold on;
for i = 1:l
  plot(x, walkers(:,i),'.');
  leg(i) = "Walker " + i;
end
hold off;
axis square
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
title('Trajactory \lambda States','FontSize', 16);
xlabel('Time (nsec)');
ylabel('\lambda');
legend(leg, 'Location', 'SouthWest');
