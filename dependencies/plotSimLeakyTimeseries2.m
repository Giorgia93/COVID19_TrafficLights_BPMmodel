function [] = plotSimLeakyTimeseries2(AL_history, daily_cases, daily_hosp, n_deaths, hosp_beds, TTIQeff_hist, scenario, tol, border)

t = 1:length(daily_cases);
% scenario = "A";
% tol = "high";
% border = "high";

nrows = size(AL_history, 1);
ncols = size(AL_history, 2);

G_props = sum(AL_history == 1)/nrows;
Y_props = sum(AL_history == 2)/nrows;
R_props = sum(AL_history == 3)/nrows;


colours = ["#66FF66" "#FFFF66" "#FF6666" "#6666FF"];
col_95prct = [1 1 1];
col_50prct = [0.87 0.87 0.87];

f = figure;
f.Position = [100 100 1200 250];
tl = tiledlayout(1,4);
% title(tl,append('Scenario ', scenario, ', ', tol, ' tolerance, ', border, ' border'))
title(tl,append(tol, ' tolerance'))

%%%% Plot daily reported cases
nexttile
prct95 = prctile(daily_cases(:,1:end-1), [2.5 97.5]);
prct50 = prctile(daily_cases(:,1:end-1), [25 75]);

hold on
for ti = 1:ncols
    area([ti,ti+1],[1, 1] .* max(prct95(2,:))*1.15, 'FaceColor', colours(4), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti) + Y_props(ti) + R_props(ti), G_props(ti) + Y_props(ti) + R_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(3), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti) + Y_props(ti), G_props(ti) + Y_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(2), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti), G_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(1), 'LineStyle', 'none')
end
fill([t(:,1:end-1), fliplr(t(:,1:end-1))], [prct95(1, :), fliplr(prct95(2, :))], col_95prct, 'FaceAlpha', 0.3, 'LineStyle', 'none')
fill([t(:,1:end-1), fliplr(t(:,1:end-1))], [prct50(1, :), fliplr(prct50(2, :))], col_95prct, 'FaceAlpha', 0.5, 'LineStyle', 'none')
plot(t(:,1:end-1), mean(daily_cases(:,1:end-1)), 'k-')
hold off
datetick('x', 'dd mmm')
ylabel("daily reported cases")
ylim([0 max(prct95(2,:))*1.1]);


%%%% Plot daily hospitalisations
nexttile
prct95 = prctile(daily_hosp, [2.5 97.5]);
prct50 = prctile(daily_hosp, [25 75]);
hold on
for ti = 1:ncols
    area([ti,ti+1],[1, 1] .* max(prct95(2,:))*1.15, 'FaceColor', colours(4), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti) + Y_props(ti) + R_props(ti), G_props(ti) + Y_props(ti) + R_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(3), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti) + Y_props(ti), G_props(ti) + Y_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(2), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti), G_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(1), 'LineStyle', 'none')
end
fill([t, fliplr(t)], [prct95(1, :), fliplr(prct95(2, :))], col_95prct, 'FaceAlpha', 0.3, 'LineStyle', 'none')
fill([t, fliplr(t)], [prct50(1, :), fliplr(prct50(2, :))], col_95prct, 'FaceAlpha', 0.5, 'LineStyle', 'none')
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
for ti = 1:ncols
    area([ti,ti+1],[1, 1] .* max(prct95(2,:))*1.15, 'FaceColor', colours(4), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti) + Y_props(ti) + R_props(ti), G_props(ti) + Y_props(ti) + R_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(3), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti) + Y_props(ti), G_props(ti) + Y_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(2), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti), G_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(1), 'LineStyle', 'none')
end
fill([t, fliplr(t)], [prct95(1, :), fliplr(prct95(2, :))], col_95prct, 'FaceAlpha', 0.3, 'LineStyle', 'none')
fill([t, fliplr(t)], [prct50(1, :), fliplr(prct50(2, :))], col_50prct, 'FaceAlpha', 0.5, 'LineStyle', 'none')
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
for ti = 1:ncols
    area([ti,ti+1],[1, 1] .* max(prct95(2,:))*1.15, 'FaceColor', colours(4), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti) + Y_props(ti) + R_props(ti), G_props(ti) + Y_props(ti) + R_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(3), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti) + Y_props(ti), G_props(ti) + Y_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(2), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti), G_props(ti)] .* max(prct95(2,:))*1.15, 'FaceColor', colours(1), 'LineStyle', 'none')
end
fill([t(:,1:end-1), fliplr(t(:,1:end-1))], [prct95(1, :), fliplr(prct95(2, :))], col_95prct, 'FaceAlpha', 0.3, 'LineStyle', 'none')
fill([t(:,1:end-1), fliplr(t(:,1:end-1))], [prct50(1, :), fliplr(prct50(2, :))], col_95prct, 'FaceAlpha', 0.5, 'LineStyle', 'none')
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