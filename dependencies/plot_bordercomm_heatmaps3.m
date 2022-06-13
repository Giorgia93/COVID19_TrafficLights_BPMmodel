% clear

% Import the data
commborder0 = readtable("results\summaries\commborderNEW.csv");
commborder0.infections = round(commborder0.infections, -3);

commborder = grpstats(commborder0,["tolerance", "iniCases","borderCases"], "mean", ...
    "DataVars", ["infections", "red", "emergency"], ...
    "VarNames",["tolerance", "iniCases", "borderCases", "count", "infections", "red", "emergency"]);
commborder.infections = round(commborder.infections, -3);
commborder.Properties.RowNames = {};

sum = grpstats(commborder0,["tolerance", "iniCases","borderCases"], "mean", ...
    "DataVars", ["infections", "cases", "hosp", "peakBeds", "deaths", ...
    "green", "yellow", "red", "emergency", "TTIQeff"], ...
    "VarNames",["tolerance", "iniCases", "borderCases", "count", ...
    "infections", "cases", "hosp", "peakBeds", "deaths", "green", "yellow", "red", "emergency", "TTIQeff"]);
sum.Properties.RowNames = {};


% TOTAL INFECTIONS

f = figure;
% title(f, "Total infections")
f.Position = [50 -100 1400 900];
% sgtitle('Total infections / year')
% nrow, ncol, [vspace, hspace], [space below, space above], [space left, space right]
[ha, ~] = tight_subplot(3, 4, [.13 .09], [.05 .05],[.05 .05]);

axes(ha(1))
% title('Total infections / year')
VLT = commborder(commborder.tolerance == "vlow", :);
h = heatmap(VLT,'iniCases','borderCases','ColorVariable','infections');
h.CellLabelFormat = '%.0f';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
title("Very low tolerance")
xlabel("initial community infections")
ylabel("Border infections / year")

axes(ha(2))
LT = commborder(commborder.tolerance == "low", :);
h = heatmap(LT,'iniCases','borderCases','ColorVariable','infections');
h.CellLabelFormat = '%.0f';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
title("Low tolerance")
xlabel("initial community infections")
ylabel("border infections / year")

axes(ha(3))
MT = commborder(commborder.tolerance == "medium", :);
h = heatmap(MT,'iniCases','borderCases','ColorVariable','infections');
h.CellLabelFormat = '%.0f';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
title("Medium tolerance")
xlabel("initial community infections")
ylabel("border infections / year")

axes(ha(4))
HT = commborder(commborder.tolerance == "HT" | commborder.tolerance == "high", :);
h = heatmap(HT,'iniCases','borderCases','ColorVariable','infections');
h.CellLabelFormat = '%.0f';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
title("High tolerance")
xlabel("initial community infections")
ylabel("border infections / year")

% PROP. TIME AT RED

% f = figure;
% title(f, "Total infections")
% f.Position = [300 300 1400 300];
% sgtitle('Proportion of time spent at Red')
% nrow, ncol, [vspace, hspace], [space below, space above], [space left, space right]
% [ha, pos] = tight_subplot(1,4,[.03 .08],[.2 .2],[.08 .08]);

axes(ha(5))
VLT = commborder(commborder.tolerance == "vlow", :);
% VLT.red = VLT.red.*100;
h = heatmap(VLT,'iniCases','borderCases','ColorVariable','red');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [ones(1, 101); (1:-0.01:0); (1:-0.01:0)]';
h.ColorLimits = [0 100];
title("Very low tolerance")
xlabel("initial community infections")
ylabel("border infections / year")

axes(ha(6))
LT = commborder(commborder.tolerance == "low", :);
% LT.red = LT.red.*100;
h = heatmap(LT,'iniCases','borderCases','ColorVariable','red');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [ones(1, 101); (1:-0.01:0); (1:-0.01:0)]';
h.ColorLimits = [0 100];
title("Low tolerance")
xlabel("initial community infections")
ylabel("border infections / year")

axes(ha(7))
MT = commborder(commborder.tolerance == "medium", :);
% MT.red = MT.red.*100;
h = heatmap(MT,'iniCases','borderCases','ColorVariable','red');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [ones(1, 101); (1:-0.01:0); (1:-0.01:0)]';
h.ColorLimits = [0 100];
title("Medium tolerance")
xlabel("initial community infections")
ylabel("border infections / year")

axes(ha(8))
HT = commborder(commborder.tolerance == "high", :);
% HT.red = HT.red.*100;
h = heatmap(HT,'iniCases','borderCases','ColorVariable','red');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [ones(1, 101); (1:-0.01:0); (1:-0.01:0)]';
h.ColorLimits = [0 100];
title("High tolerance")
xlabel("initial community infections")
ylabel("border infections / year")


% PROP TIME AT EMERGENCY

% f = figure;
% % title(f, "Total infections")
% f.Position = [300 300 1400 300];
% sgtitle('Proportion of time spent at Emergency')
% nrow, ncol, [vspace, hspace], [space below, space above], [space left, space right]
% [ha, pos] = tight_subplot(1,4,[.03 .08],[.2 .2],[.08 .08]);

axes(ha(9))
VLT = commborder(commborder.tolerance == "vlow", :);
% VLT.emergency = VLT.emergency .*100;
h = heatmap(VLT,'iniCases','borderCases','ColorVariable','emergency');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [linspace(1, 0.5, 100); linspace(1, 0, 100); linspace(1, 1, 100)]';
h.ColorLimits = [0 100];
title("Very low tolerance")
xlabel("initial community infections")
ylabel("border infections / year")

axes(ha(10))
LT = commborder(commborder.tolerance == "low", :);
% LT.emergency = LT.emergency .*100;
h = heatmap(LT,'iniCases','borderCases','ColorVariable','emergency');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [linspace(1, 0.5, 100); linspace(1, 0, 100); linspace(1, 1, 100)]';
h.ColorLimits = [0 100];
title("Low tolerance")
xlabel("Proportion of time at Emergency")
ylabel("border infections / year")

axes(ha(11))
MT = commborder(commborder.tolerance == "medium", :);
% MT.emergency = MT.emergency .*100;
h = heatmap(MT,'iniCases','borderCases','ColorVariable','emergency');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [linspace(1, 0.5, 100); linspace(1, 0, 100); linspace(1, 1, 100)]';
h.ColorLimits = [0 100];
title("Medium tolerance")
xlabel("initial community infections")
ylabel("border infections / year")

axes(ha(12))
HT = commborder(commborder.tolerance == "high", :);
% HT.emergency = HT.emergency .*100;
h = heatmap(HT,'iniCases','borderCases','ColorVariable','emergency');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.borderCases = categorical(h.SourceTable.borderCases);
h.SourceTable.borderCases = reordercats(h.SourceTable.borderCases,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [linspace(1, 0.5, 100); linspace(1, 0, 100); linspace(1, 1, 100)]';
h.ColorLimits = [0 100];
title("High tolerance")
xlabel("initial community infections")
ylabel("border infections / year")
