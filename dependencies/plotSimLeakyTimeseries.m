function [] = plotSimLeakyTimeseries(AL_history, daily_cases, daily_hosp, n_deaths, hosp_beds, TTIQeff_hist, scenario, tol, border)

t = 1:length(daily_cases);
% scenario = "A";
% tol = "high";
% border = "high";

AlHist = round(mean(AL_history));

colours = ["#66FF66" "#FFFF66" "#FF6666" "#6666FF"];
col_95prct = [1 1 1];
col_50prct = [0.87 0.87 0.87];

tl = tiledlayout(2,2);
title(tl,append('Scenario ', scenario, ', ', tol, ' tolerance, ', border, ' border'))

%%%% Plot daily reported cases
nexttile
prct95 = prctile(daily_cases(:,1:end-1), [2.5 97.5]);
prct50 = prctile(daily_cases(:,1:end-1), [25 75]);
hold on
for ti = t(1:end - 1)
    area([ti,ti+1],[max(prct95(2,:))*1.1, max(prct95(2,:))*1.1], 'FaceColor', colours(AlHist(ti)), 'LineStyle', 'none')
end
fill([t(:,1:end-1), fliplr(t(:,1:end-1))], [prct95(1, :), fliplr(prct95(2, :))], col_95prct, 'FaceAlpha', 0.8, 'LineStyle', 'none')
fill([t(:,1:end-1), fliplr(t(:,1:end-1))], [prct50(1, :), fliplr(prct50(2, :))], col_95prct, 'FaceAlpha', 1, 'LineStyle', 'none')
plot(t(:,1:end-1), mean(daily_cases(:,1:end-1)), 'k')
hold off
datetick('x', 'dd mmm')
ylabel("daily reported cases")
ylim([0 max(prct95(2,:))*1.1]);

%%%% Plot daily hospitalisations
nexttile
prct95 = prctile(daily_hosp, [2.5 97.5]);
prct50 = prctile(daily_hosp, [25 75]);
hold on
for ti = t(1:end - 1)
    area([ti,ti+1],[max(prct95(2,:))*1.1, max(prct95(2,:))*1.1], 'FaceColor', colours(AlHist(ti)), 'LineStyle', 'none')
end
fill([t, fliplr(t)], [prct95(1, :), fliplr(prct95(2, :))], col_95prct, 'LineStyle', 'none')
fill([t, fliplr(t)], [prct50(1, :), fliplr(prct50(2, :))], col_50prct, 'LineStyle', 'none')
plot(t, mean(daily_hosp), 'k')
hold off
datetick('x', 'dd mmm')
ylabel("new daily hospitalisations")
ylim([0 max(prct95(2,:))*1.1]);

%%%% Plot cumulative deaths
cDeaths = zeros(size(n_deaths));
for i = 1:length(t); cDeaths(:, i) = sum(n_deaths(:, 1:i), 2); end
nexttile
prct95 = prctile(cDeaths, [2.5 97.5]);
prct50 = prctile(cDeaths, [25 75]);
hold on
for ti = t(1:end - 1)
    area([ti,ti+1],[max(prct95(2,:))*1.1, max(prct95(2,:))*1.1], 'FaceColor', colours(AlHist(ti)), 'LineStyle', 'none')
end
fill([t, fliplr(t)], [prct95(1, :), fliplr(prct95(2, :))], col_95prct, 'LineStyle', 'none')
fill([t, fliplr(t)], [prct50(1, :), fliplr(prct50(2, :))], col_50prct, 'LineStyle', 'none')
plot(t, mean(cDeaths), 'k')
hold off
datetick('x', 'dd mmm')
ylabel("cumulative deaths")
ylim([0 max(prct95(2,:))*1.1]);

%%%% Plot hospital beds occupied
nexttile
prct95 = prctile(hosp_beds(:,1:end-1), [2.5 97.5]);
prct50 = prctile(hosp_beds(:,1:end-1), [25 75]);
hold on
for ti = t(1:end - 1)
    area([ti,ti+1],[max(prct95(2,:))*1.1, max(prct95(2,:))*1.1], 'FaceColor', colours(AlHist(ti)), 'LineStyle', 'none')
end
fill([t(:,1:end-1), fliplr(t(:,1:end-1))], [prct95(1, :), fliplr(prct95(2, :))], col_95prct, 'LineStyle', 'none')
fill([t(:,1:end-1), fliplr(t(:,1:end-1))], [prct50(1, :), fliplr(prct50(2, :))], col_50prct, 'LineStyle', 'none')
plot(t(:,1:end-1), mean(hosp_beds(:,1:end-1)), 'k')
hold off
datetick('x', 'dd mmm')
ylabel("hospital beds occupied")
ylim([0 max(prct95(2,:))*1.1]);


% relTransAL = [0.9, 0.8, 0.7, 0.4];
% Rv = 1.7473;
% Reff_history = (1 - mean(TTIQeff_hist)) .* relTransAL(round(mean(AL_history))) .* Rv;
% figure
% hold on
% % fill([t(:,1:end-1), fliplr(t(:,1:end-1))], [prct95(1, :), fliplr(prct95(2, :))], col_95prct, 'LineStyle', 'none')
% % fill([t(:,1:end-1), fliplr(t(:,1:end-1))], [prct50(1, :), fliplr(prct50(2, :))], col_50prct, 'LineStyle', 'none')
% plot(t(:,1:end-1), Reff_history(1:end-1))
% hold off
% datetick('x', 'dd mmm')
% ylabel("Reff")


end