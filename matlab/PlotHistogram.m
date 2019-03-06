clear all;
close all;
figure(1);
ffxHist = load('histogram.txt');
% plot(ffxHist(:,2),ffxHist(:,6),'--r',ffxHist(:,2),ffxHist(:,7),'-b',ffxHist(:,2),ffxHist(:,8),ffxHist(:,2),ffxHist(:,9),ffxHist(:,2),ffxHist(:,10),ffxHist(:,2),ffxHist(:,11))
% legend('<dU/dL>', 'g(L)', 'f(L,<dU/dL>)', 'Bias', 'Free Energy', 'Free Energy + Bias')
plot(ffxHist(:,2),ffxHist(:,6),'-k',ffxHist(:,2),ffxHist(:,10),'--g',ffxHist(:,2),ffxHist(:,11),'.b');
legend('<dU/dL>', 'Free Energy', 'Free Energy + Bias', 'FontSize', 14, 'FontName', 'Times New Roman');
title('Sampling Bias and Free Energy','FontSize', 16, 'FontName', 'Times New Roman');
xlabel('Lambda','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Energy (kcal/mol)','FontSize', 14, 'FontName', 'Times New Roman');
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times New Roman';

figure(2)
ffxPMF = load('pmf.txt');
y = ffxPMF(:,1);
x = ffxHist(:,2);
data = ffxPMF(:,2:202);
contour(x,y,data,20);
colorbar();
hold on;
plot(ffxHist(:,2),ffxHist(:,6),'-k');
title('Potential of Mean Force');
xlabel('Lambda');
ylabel('dU/dL (kcal/mol)');
legend('Total Bias PMF','<dU/dL>');
hold off;

figure(3)
ffxBias = load('pmf.2D.txt');
y = ffxBias(:,1);
data = ffxBias(:,2:202);
contour(x,y,data,5);
colorbar();
hold on;
plot(ffxHist(:,2),ffxHist(:,6),'-k');
title('2D Bias Potential of Mean Force');
xlabel('Lambda');
ylabel('dU/dL (kcal/mol)');
legend('2D Bias PMF','<dU/dL>');
hold off;