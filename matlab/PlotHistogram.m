clear all;
close all;
figure(1);
ffxHist = load('histogram.txt');
plot(ffxHist(:,2),ffxHist(:,5),'-k',ffxHist(:,2),ffxHist(:,9),'--g',ffxHist(:,2),ffxHist(:,10),'.b');
axis square
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
title('Sampling Bias and Free Energy','FontSize', 16);
xlabel('Lambda');
ylabel('Energy (kcal/mol)');
legend('<dU/dL>', 'Free Energy', 'Free Energy + Bias');

figure(2);
ffxPMF = load('pmf.txt');
y = ffxPMF(:,1);
x = ffxHist(:,2);
data = ffxPMF(:,2:202);
contour(x,y,data,20);
colorbar();
hold on;
plot(ffxHist(:,2),ffxHist(:,5),'-k');
axis square
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
title('Potential of Mean Force','FontSize', 16);
xlabel('Lambda');
ylabel('dU/dL (kcal/mol)');
legend('Total Bias PMF','<dU/dL>');
hold off;

figure(3);
ffxBias = load('pmf.2D.txt');
y = ffxBias(:,1);
data = ffxBias(:,2:202);
contour(x,y,data,5);
colorbar();
hold on;
plot(ffxHist(:,2),ffxHist(:,5),'-k');
axis square
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
title('2D Bias Potential of Mean Force','FontSize', 16);
xlabel('Lambda');
ylabel('dU/dL (kcal/mol)');
legend('2D Bias PMF','<dU/dL>');
hold off;