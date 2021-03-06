load spreadsheetdata.mat

distmeters = distfeet * 0.3;

%pmf = flipud(rssipmf); % 61x20 for 61 attemiations bins
pmf = rssipmf;
figure;
imagesc(pmf);
xt = xticks;
xticklabels(distmeters(xt));
xlabel('distance (meters)')
yt = yticks;
yticklabels(dBms(yt));
ylabel('RSSI')
set(gca,'YDir','normal')
print('rss_dist_pmf_heatmap', '-dpng')

Y = rssivals'; % 240 x 20 for distance=1:20 feet
figure;
boxplot(Y)
x = distmeters;
xticklabels(x)
xlabel('distance (meters)')
ylabel('RSSI')
print('rss_dist_boxplot', '-dpng')


% Tom Lovett paper, Guassian noise model, p20
% https://arxiv.org/abs/2007.05057
x = distmeters;
rss_pred = -8.69*log(x) - 67.9;
sigma = sqrt(97.03);
err = sigma*ones(size(x));
figure;
errorbar(x, rss_pred, err)
xlabel('distance (meters)')
ylabel('predicted RSSI')
title('Lovett')
print('rss_dist_pred_lovett', '-dpng')

% Fit ridge regression model
Y = rssivals'; 
y = Y(:);
N = size(Y,1);
X = repmat(distmeters, N, 1);
x = X(:);

indices = find(isnan(y) == 0);
[I,J] = ind2sub(size(y),indices);
ndx = I;
y = y(ndx);
x = x(ndx);

w = ridge(y, log(x), 1e-5, 0) % -63.00, -5.85
rss_pred = w(1) + w(2)*log(x);
N = length(y);
mse = sum((y - rss_pred).^2)/N
sigma = sqrt(mse);
err = sigma*ones(size(x));
figure;
errorbar(x, rss_pred, err)
xlabel('distance (meters)')
ylabel('predicted RSSI')
title('fitted')
print('rss_dist_pred_ridge', '-dpng')



Y = rssivals'; % 240 x 20 for distance=1:20 feet
figure; hold on;
boxplot(Y)
n = size(Y,2);
x = distmeters;
rss_pred = w(1) + w(2)*log(x);
%rss_pred = -8.69*log(x) - 67.9;
plot(1:n, rss_pred);
xticklabels(distmeters)
xlabel('distance (meters)')
ylabel('RSSI')
print('rss_dist_boxplot_pred', '-dpng')
