figure

cats = categorical({'infections', 'cases', 'hospitalisations', 'peakbeds', 'deaths', 'prop. time at R', 'prop. time at E'});
cats = reordercats(cats,{'infections', 'cases', 'hospitalisations', 'peakbeds', 'deaths', 'prop. time at R', 'prop. time at E'});
lowtol = [0.751775958	0.794670484	2.139980697	1.669319167	0.863429186			-0.625113464	3.538344434];
medtol = [1.076077389	1.162650639	2.697531901	1.833172961	1.262701817			-0.254386306	7.681609195];
hightol = [0.468565607	0.528830587	1.609114772	1.242134081	0.562499612			-0.005028902	7.09928401];
y = [lowtol', medtol', hightol'];


y = [-1.95 	-1.13 	-3.14 	-3.33 	-2.73 	-0.33	-0.15;
-0.44 	-0.45 	-0.77 	-1.76 	-0.72 	-0.28	-0.05;
-0.36 	-0.36 	-0.67 	-1.24 	-0.58 	-0.20	-0.06
]';



b = barh(cats, y,'BaseValue',0);
set(gca, 'YDir','reverse')
yticklabels(cats)
% leg = legend({'Low tolerance', 'Medium tolerance','High tolerance'}, 'location', 'northeastoutside');
% title(leg,'Scenario')
xlim([-4, 0])
% title({'Percentage impact with', 'lowered vaccination efficacy'});
set(gca,'xtick',[])
set(gca,'xticklabel',[])

for i = 1:3
    xtips = b(i).YEndPoints + (0.002 .* b(i).YEndPoints./abs(b(i).YEndPoints));
    xtips(xtips < 0) = xtips(xtips < 0) - 0.5;
    ytips = b(i).XEndPoints;
    pluses = strings(1, length(xtips));
    pluses(xtips > 0) = "+";
    labels = append(pluses, string(round(100 .* (b(i).YData))), '%');
    
    text(xtips,ytips,labels,'VerticalAlignment','middle')
end

