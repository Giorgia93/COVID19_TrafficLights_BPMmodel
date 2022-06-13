function [] = plotPropAL()

scenariofolder = "A-2";
scenario = "A";
tol = "high";
border = "high";
date = "20211118";
AL_history = table2array(readtable(append('results/timeseries/scenario', scenariofolder, '/ALhistory_scenario', ...
    scenario, '_', tol, 'tol_', border, 'border_', date, '.csv')));

nrows = size(AL_history, 1);
ncols = size(AL_history, 2);

G_props = sum(AL_history == 1)/nrows;
Y_props = sum(AL_history == 2)/nrows;
R_props = sum(AL_history == 3)/nrows;

colours = ["#66FF66" "#FFFF66" "#FF6666" "#6666FF"];

figure
hold on
for ti = 1:ncols
    area([ti,ti+1],[1, 1], 'FaceColor', colours(4), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti) + Y_props(ti) + R_props(ti), G_props(ti) + Y_props(ti) + R_props(ti)], 'FaceColor', colours(3), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti) + Y_props(ti), G_props(ti) + Y_props(ti)], 'FaceColor', colours(2), 'LineStyle', 'none')
    area([ti,ti+1],[G_props(ti), G_props(ti)], 'FaceColor', colours(1), 'LineStyle', 'none')
end
hold off
xlim([0,ncols])
datetick('x', 'dd mmm')
ylabel("proportion of simulations")

end