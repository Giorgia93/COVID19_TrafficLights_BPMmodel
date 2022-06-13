function [] = plotSimLeaky(t, nActiveDetected, nHosp, nHospBeds, nDeaths)

tiledlayout(2,2);

%%%% Plot daily reported cases
nexttile
plot(t, sum(nActiveDetected, 1))
datetick('x', 'dd mmm')
ylabel("daily reported cases")


%%%% Plot daily hospitalisations
nexttile
plot(t, sum(nHosp, 2))
datetick('x', 'dd mmm')
ylabel("new daily hospitalisations")

%%%% Plot cumulative deaths
cDeaths = zeros(length(t), 16);
nDeaths2 = sum(nDeaths, 2);
for i = 1:length(t); cDeaths(i, :) = sum(nDeaths2(1:i, :)); end
nexttile
plot(t, cDeaths)
datetick('x', 'dd mmm')
ylabel("cumulative deaths")

%%%% Plot hospital beds occupied
nexttile
plot(t, nHospBeds)
datetick('x', 'dd mmm')
ylabel("hospital beds occupied")


end