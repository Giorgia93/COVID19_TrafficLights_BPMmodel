%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 15);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["scen", "tol", "Comm0", "Border", "infections", "cases", "hosp", "peakBeds", "deaths", "G", "Y", "R", "E", "TTIQeff", "sims"];
opts.VariableTypes = ["categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["scen", "tol"], "EmptyFieldRule", "auto");

% Import the data
commborder = readtable("\\file\Usersg$\gva31\Home\My Documents\GiorgiasSandpit\GiorgiasSimLeaky\results\summaries\commborder2.csv", opts);


%% Clear temporary variables
clear opts
% close all


%% TOTAL INFECTIONS

f = figure;
% title(f, "Total infections")
f.Position = [300 300 1400 300];
sgtitle('Total infections / year')
% nrow, ncol, [vspace, hspace], [space below, space above], [space left, space right]
[ha, pos] = tight_subplot(1,4,[.03 .08],[.2 .2],[.08 .08]);

axes(ha(1))
VLT = commborder(commborder.tol == "VLT" | commborder.tol == "vlow", :);
h = heatmap(VLT,'Comm0','Border','ColorVariable','infections');
h.CellLabelFormat = '%i';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
title("Very low tolerance")
xlabel("Initial community infections")
ylabel("Border infections / year")

axes(ha(2))
LT = commborder(commborder.tol == "LT" | commborder.tol == "low", :);
h = heatmap(LT,'Comm0','Border','ColorVariable','infections');
h.CellLabelFormat = '%i';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
title("Low tolerance")
xlabel("Initial community infections")
ylabel("Border infections / year")

axes(ha(3))
MT = commborder(commborder.tol == "MT" | commborder.tol == "medium", :);
h = heatmap(MT,'Comm0','Border','ColorVariable','infections');
h.CellLabelFormat = '%i';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
title("Medium tolerance")
xlabel("Initial community infections")
ylabel("Border infections / year")

axes(ha(4))
HT = commborder(commborder.tol == "HT" | commborder.tol == "high", :);
h = heatmap(HT,'Comm0','Border','ColorVariable','infections');
h.CellLabelFormat = '%i';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
title("High tolerance")
xlabel("Initial community infections")
ylabel("Border infections / year")

%% PROP. TIME AT RED

f = figure;
% title(f, "Total infections")
f.Position = [300 300 1400 300];
sgtitle('Proportion of time spent at Red')
% nrow, ncol, [vspace, hspace], [space below, space above], [space left, space right]
[ha, pos] = tight_subplot(1,4,[.03 .08],[.2 .2],[.08 .08]);

axes(ha(1))
VLT = commborder(commborder.tol == "VLT" | commborder.tol == "vlow", :);
VLT.R = VLT.R.*100;
h = heatmap(LT,'Comm0','Border','ColorVariable','R');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [ones(1, 101); (1:-0.01:0); (1:-0.01:0)]';
h.ColorLimits = [0 100];
title("Very low tolerance")
xlabel("Initial community infections")
ylabel("Border infections / year")

axes(ha(2))
LT = commborder(commborder.tol == "LT" | commborder.tol == "low", :);
LT.R = LT.R.*100;
h = heatmap(LT,'Comm0','Border','ColorVariable','R');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [ones(1, 101); (1:-0.01:0); (1:-0.01:0)]';
h.ColorLimits = [0 100];
title("Low tolerance")
xlabel("Initial community infections")
ylabel("Border infections / year")

axes(ha(3))
MT = commborder(commborder.tol == "MT" | commborder.tol == "medium", :);
MT.R = MT.R.*100;
h = heatmap(MT,'Comm0','Border','ColorVariable','R');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [ones(1, 101); (1:-0.01:0); (1:-0.01:0)]';
h.ColorLimits = [0 100];
title("Medium tolerance")
xlabel("Initial community infections")
ylabel("Border infections / year")

axes(ha(4))
HT = commborder(commborder.tol == "HT" | commborder.tol == "high", :);
HT.R = HT.R.*100;
h = heatmap(HT,'Comm0','Border','ColorVariable','R');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [ones(1, 101); (1:-0.01:0); (1:-0.01:0)]';
h.ColorLimits = [0 100];
title("High tolerance")
xlabel("Initial community infections")
ylabel("Border infections / year")


%% PROP TIME AT EMERGENCY

f = figure;
% title(f, "Total infections")
f.Position = [300 300 1400 300];
sgtitle('Proportion of time spent at Emergency')
% nrow, ncol, [vspace, hspace], [space below, space above], [space left, space right]
[ha, pos] = tight_subplot(1,4,[.03 .08],[.2 .2],[.08 .08]);

axes(ha(1))
VLT = commborder(commborder.tol == "VLT" | commborder.tol == "vlow", :);
VLT.E = VLT.E .*100;
h = heatmap(VLT,'Comm0','Border','ColorVariable','E');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [linspace(1, 0.5, 100); linspace(1, 0, 100); linspace(1, 1, 100)]';
h.ColorLimits = [0 100];
title("Very low tolerance")
xlabel("Initial community infections")
ylabel("Border infections / year")

axes(ha(2))
LT = commborder(commborder.tol == "LT" | commborder.tol == "low", :);
LT.E = LT.E .*100;
h = heatmap(LT,'Comm0','Border','ColorVariable','E');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [linspace(1, 0.5, 100); linspace(1, 0, 100); linspace(1, 1, 100)]';
h.ColorLimits = [0 100];
title("Low tolerance")
xlabel("Proportion of time at Emergency")
ylabel("Border infections / year")

axes(ha(3))
MT = commborder(commborder.tol == "MT" | commborder.tol == "medium", :);
MT.E = MT.E .*100;
h = heatmap(MT,'Comm0','Border','ColorVariable','E');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [linspace(1, 0.5, 100); linspace(1, 0, 100); linspace(1, 1, 100)]';
h.ColorLimits = [0 100];
title("Medium tolerance")
xlabel("Initial community infections")
ylabel("Border infections / year")

axes(ha(4))
HT = commborder(commborder.tol == "HT" | commborder.tol == "high", :);
HT.E = HT.E .*100;
h = heatmap(HT,'Comm0','Border','ColorVariable','E');
h.CellLabelFormat = '%.0f%%';
h.SourceTable.Border = categorical(h.SourceTable.Border);
h.SourceTable.Border = reordercats(h.SourceTable.Border,{'20000', '10000', '5000', '500', '0'});
h.Colormap = [linspace(1, 0.5, 100); linspace(1, 0, 100); linspace(1, 1, 100)]';
h.ColorLimits = [0 100];
title("High tolerance")
xlabel("Initial community infections")
ylabel("Border infections / year")
