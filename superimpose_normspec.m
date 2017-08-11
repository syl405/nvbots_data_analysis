figure
%% narrow allowable range
subplot(1,2,1)
hold on

mu = 0;              % mean
sigma = 0.2;          % std
x = linspace(mu-5*sigma, mu+5*sigma, 500);
low2 = -0.2;
high2 = 0.2;
y = normpdf(x, mu, sigma);
xci2 = [linspace(mu-5*sigma, low2); linspace(high2, mu+5*sigma)];
yci2 = normpdf(xci2, mu, sigma);
precomp_fh = plot(x, y, '-r', 'LineWidth', 1.5);
%patch(x, y, [1 1 1],'facealpha',0.5)
patch([xci2(1,:) low2], [yci2(1,:) 0], [1 0 0],'facealpha',0.5)
patch([high2 xci2(2,:)], [0 yci2(2,:)], [1 0 0],'facealpha',0.5)

mu = 0;              % mean
sigma = 0.3582;          % std
x = linspace(mu-5*sigma, mu+5*sigma, 500);
low2 = -0.2;
high2 = 0.2;
y = normpdf(x, mu, sigma);
xci2 = [linspace(mu-5*sigma, low2); linspace(high2, mu+5*sigma)];
yci2 = normpdf(xci2, mu, sigma);
postcomp_fh = plot(x, y, '-g', 'LineWidth', 1.5);
%patch(x, y, [1 1 1],'facealpha',0.5)
patch([xci2(1,:) low2], [yci2(1,:) 0], [0 1 0],'facealpha',0.5)
patch([high2 xci2(2,:)], [0 yci2(2,:)], [0 1 0],'facealpha',0.5)
xlabel('Nozzle offset (mm)')
ylabel('Probability Density')

legend([precomp_fh postcomp_fh],{'Before Compensation' 'After Compensation'})
title('Allowable nozzle offset range = \pm 0.2 mm')

%% wide allowable range
subplot(1,2,2)
hold on

mu = 0;              % mean
sigma = 0.2;          % std
x = linspace(mu-5*sigma, mu+5*sigma, 500);
low1 = -0.5;
high1 = 0.5;
y = normpdf(x, mu, sigma);
xci1 = [linspace(mu-5*sigma, low1); linspace(high1, mu+5*sigma)];
yci1 = normpdf(xci1, mu, sigma);
precomp_fh = plot(x, y, '-r', 'LineWidth', 1.5);
%patch(x, y, [1 1 1],'facealpha',0.5)
patch([xci1(1,:) low1], [yci1(1,:) 0], [1 0 0],'facealpha',0.5)
patch([high1 xci1(2,:)], [0 yci1(2,:)], [1 0 0],'facealpha',0.5)

mu = 0;              % mean
sigma = 0.3582;          % std
x = linspace(mu-5*sigma, mu+5*sigma, 500);
low1 = -0.5;
high1 = 0.5;
y = normpdf(x, mu, sigma);
xci1 = [linspace(mu-5*sigma, low1); linspace(high1, mu+5*sigma)];
yci1 = normpdf(xci1, mu, sigma);
postcomp_fh = plot(x, y, '-g', 'LineWidth', 1.5);
%patch(x, y, [1 1 1],'facealpha',0.5)
patch([xci1(1,:) low1], [yci1(1,:) 0], [0 1 0],'facealpha',0.5)
patch([high1 xci1(2,:)], [0 yci1(2,:)], [0 1 0],'facealpha',0.5)
xlabel('Nozzle offset (mm)')
ylabel('Probability Density')

legend([precomp_fh postcomp_fh],{'Before Compensation' 'After Compensation'})
title('Allowable nozzle offset range = \pm 0.5 mm')