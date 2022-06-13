%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   BRANCHING PROCESS MODEL
%                   COVID RESPONSE STRATEGY
%                    TRAFFIC LIGHTS SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc


%------------------------- SIMULATION PARAMETERS --------------------------

nsims = 2; % Number of simulations for each combination of params
today = datestr(now,'yyyymmdd');

cov_scenarios = "90% of 12+"; %["90% of 12+", "90% of 5+", "95% of 12+", "95% of 5+"];
cov = 0.90; %[0.9, 0.9, 0.95, 0.95];
covage = 12; %[12, 5, 12, 5];

tol_levels = ["low", "medium", "high"]; %["low", "medium", "high"];

% Testing probability of clinical cases. Baseline = 0.45
pTest_clinical = [0.45]; % [0.45, 0.35]

% Number of initial community infections. Baseline = 2000
infections0 = [2000];%[0, 500, 1000, 2000, 10000]; %[2000, 10000];
iAL0 = 3; % initial alert level

% Number of border cases throughout the year. Baseline = 5000
border_seeds = [5000];%[0, 500, 5000, 10000, 20000]; %[500, 5000]; [10000, 20000]

% Relative transmission at each alert level. Baseline = 'high' [0.9, 0.8, 0.7, 0.4]
reltrans_AL_levels = ["high"]; %["high"]; %["high", "low"]
reltrans_AL = [[0.9, 0.8, 0.7, 0.4]]; %[[0.9, 0.8, 0.7, 0.4]]; %

% Vaccine efficacy against transmission, disease and death. Baseline = 'high' [0.7, 0.5, 0.8]
VE2_levels = ["high"]; % ["high", "low"]
VE2 = [[0.7, 0.5, 0.8]]; % [[0.7, 0.5, 0.8]; [0.5, 0.4, 0.8]]

% Contact tracing capacity. Baseline = 100
trace_capacity = [100]; % [100, 250]

% Probability of detecting through contact tracing. Baseline = 0.7
p_trace = [0.7]; %[0.7, 0.3]

% Border cases infectiousness Baseline = 0.2
border_inf = [0.2]; %[0.7, 0.3]

% Isolated comm. cases infectiousness Baseline = 0
isol_inf = [0]; %[0.7, 0.3]

% If true, hospitalisation rates are changed to OR=2. Baseline = false
HR_bool = true;


scenarios_combination = combvec(1:size(tol_levels, 2),...
    pTest_clinical, infections0, border_seeds, ...
    1:size(reltrans_AL_levels, 2), 1:size(VE2_levels, 1), trace_capacity,...
    p_trace, border_inf, isol_inf);
[nvar, nscenarios] = size(scenarios_combination);



plot_timeseries = false;


% Load contact data and pop dist data from files
[C, tmp, popSizeData] = loadData('nzcontmatrix.xlsx', 'nzpopdist.xlsx', 'nzpopdist.xlsx');

% Reject simulations if too far from current case data (NOT CURRENTLY USING)
earlyReject.tData = [];
earlyReject.nCasesData = [];
earlyReject.threshold = inf;


%-------------------------- SIMULATION LOOPS ------------------------------
start_time = datestr(now, 'hh:mm');


% ------ Initialise results tables ------
daily_inf = cell(nscenarios, 1);
daily_cases = cell(nscenarios, 1);
daily_hosp = cell(nscenarios, 1);
n_deaths = cell(nscenarios, 1);
hosp_beds = cell(nscenarios, 1);
AL_history = cell(nscenarios, 1);
TTIQeff_history = cell(nscenarios, 1);
ssummary = cell(nscenarios, 1);

parfor jj = 1:nscenarios
    
    scenario_letter = char(jj + 64);
    
    itol_level = scenarios_combination(1, jj); % for each tol. level
    iptest_clin = scenarios_combination(2, jj); % for each pTestClin level
    icomm_seeds = scenarios_combination(3, jj); %initial community infections
    iborder_cases = scenarios_combination(4, jj); % border level
    ireltrans_AL = reltrans_AL(scenarios_combination(5, jj), :);
    iVE2 = VE2(scenarios_combination(6, jj), :);
    itrace_capacity = scenarios_combination(7, jj);
    ip_trace = scenarios_combination(8, jj);
    iborder_inf = scenarios_combination(9, jj);
    iisol_inf = scenarios_combination(10, jj);
    
    imin_vax_age = covage; %
    icov = cov; %coverage
    
    % ------ Initialise results tables ------
    daily_inf_temp = zeros(nsims, 361);
    daily_cases_temp = zeros(nsims, 361);
    daily_hosp_temp = zeros(nsims, 361);
    n_deaths_temp = zeros(nsims, 361);
    hosp_beds_temp = zeros(nsims, 361);
    AL_history_temp = zeros(nsims, 361);
    TTIQeff_history_temp = zeros(nsims, 361);
    
    scenarioID = strings(nsims, 1);
    cov2 = zeros(nsims, 1);
    minVaxAge = zeros(nsims, 1);
    tolerance = strings(nsims, 1);
    iniCases = zeros(nsims, 1);
    borderCases = zeros(nsims, 1);
    pTestClin = zeros(nsims, 1);
    controlEfficacy = strings(nsims, 1);
    VE = strings(nsims, 1);
    traceCap = zeros(nsims, 1);
    pTrace = zeros(nsims, 1);
    borderInf = zeros(nsims, 1);
    simN = zeros(nsims, 1);
    infections = zeros(nsims, 1);
    cases = zeros(nsims, 1);
    hosp = zeros(nsims, 1);
    peakBeds = zeros(nsims, 1);
    deaths = zeros(nsims, 1);
    green = zeros(nsims, 1);
    yellow = zeros(nsims, 1);
    red = zeros(nsims, 1);
    emergency = zeros(nsims, 1);
    TTIQeff = zeros(nsims, 1);
    
    summary = table(scenarioID, cov2, minVaxAge, tolerance, iniCases, ...
        borderCases, pTestClin, controlEfficacy, VE, traceCap, pTrace,...
        borderInf, simN, infections, cases, hosp, ...
        peakBeds, deaths, green, yellow, red, emergency, TTIQeff);
    
    
    
    for i = 1:nsims % run sims
        fprintf("Combination %s (%i/%i), Sim %i/%i\n", scenario_letter, jj, nscenarios, i, nsims)
        
        % Get parameters
        par = getParParallel(HR_bool, iisol_inf, iborder_inf, ip_trace, itrace_capacity, iVE2, ...
            ireltrans_AL, iptest_clin, iAL0, icomm_seeds, imin_vax_age, ...
            icov, iborder_cases, itol_level, C, tmp, popSizeData);
        
        % Run simulation
        verbose = false;
        [nInfected, nCases, nHosp, nHospBeds, nDeaths, ...
            AL_timeseries, AL_props, TTIQeffv, TTIQeff_time] = runSimLeaky_LessVerbose(par, earlyReject, verbose);
        
        % Store data on selected scenarios for timeseries plots only
        %if ismember(scenario_counter, timeseriesdata_scenarios_tosave)
        daily_inf_temp(i, :) = sum(nInfected, 2);
        daily_cases_temp(i, :) = sum(nCases, 2);
        daily_hosp_temp(i, :) = sum(nHosp, 2);
        n_deaths_temp(i, :) = sum(nDeaths, 2);
        hosp_beds_temp(i, :) = nHospBeds;
        AL_history_temp(i, :) = AL_timeseries;
        TTIQeff_history_temp(i, :) = TTIQeff_time;
        %end
        
        % Fill summary results table
        summary.scenarioID(i) = scenario_letter;
        summary.cov2(i) = icov;
        summary.minVaxAge(i) = imin_vax_age;
        summary.tolerance(i) = par.tol_levels(itol_level);
        summary.iniCases(i) = icomm_seeds;
        summary.borderCases(i) = iborder_cases;
        summary.pTestClin(i) = iptest_clin;
        summary.controlEfficacy(i) = reltrans_AL_levels(scenarios_combination(5, jj));
        summary.VE(i)= VE2_levels(scenarios_combination(6, jj), :);
        summary.traceCap(i) = itrace_capacity;
        summary.pTrace(i) = ip_trace;
        summary.borderInf(i) = iborder_inf;
        summary.isolInf(i) = iisol_inf;
        
        summary.simN(i) = i;
        
        summary.infections(i) = sum(nInfected, 'all');
        summary.cases(i) = sum(nCases, 'all');
        summary.hosp(i) = sum(nHosp, 'all');
        summary.peakBeds(i) = max(nHospBeds);
        summary.deaths(i) = sum(nDeaths, 'all');
        summary.green(i) = round(100*AL_props(1));
        summary.yellow(i) = round(100*AL_props(2));
        summary.red(i) = round(100*AL_props(3));
        summary.emergency(i) = round(100*AL_props(4));
        summary.TTIQeff(i) = round(100*TTIQeffv);
        
        % Write timeseries results for current border and tol levels for plots
        %if ismember(scenario_counter, timeseriesdata_scenarios_tosave)
        writeSimLeakyTimeseries(daily_inf_temp, daily_cases_temp, daily_hosp_temp, n_deaths_temp, ...
            hosp_beds_temp, AL_history_temp, TTIQeff_history_temp, scenario_letter, ...
            par.tol_levels(itol_level), today)
        %end
    end
    
    daily_inf{jj} = daily_inf_temp;
    daily_cases{jj} = daily_cases_temp;
    daily_hosp{jj} = daily_hosp_temp;
    n_deaths{jj} = n_deaths_temp;
    hosp_beds{jj} = hosp_beds_temp;
    AL_history{jj} = AL_history_temp;
    TTIQeff_history{jj} = TTIQeff_history_temp;
    ssummary{jj} = summary;
    
end

summary_table = vertcat(ssummary{:});
filename1 = append('results/summaries/simLeaky_summary_', today, '.csv');
writetable(summary_table, filename1)



% Display start and end time of simulations
disp(append('Simulations started: ', start_time))
disp(append('Simulations ended: ', datestr(now, 'hh:mm')))


%------------------------------- PLOTS ------------------------------------
if plot_timeseries == true
    scenario = "B";
    tol = "medium";
    border = "low";
    date = today; %"20220427";
    daily_inf = table2array(readtable(append('results/timeseries/dailyinf_scenario', ...
        scenario, '_', tol, 'tol_', date, '.csv')));
    daily_cases = table2array(readtable(append('results/timeseries/dailycases_scenario', ...
        scenario, '_', tol, 'tol_', date, '.csv')));
    daily_hosp =table2array(readtable(append('results/timeseries/dailyhosp_scenario', ...
        scenario, '_', tol, 'tol_', date, '.csv')));
    cumul_deaths = table2array(readtable(append('results/timeseries/cumuldeaths_scenario', ...
        scenario, '_', tol, 'tol_', date, '.csv')));
    hosp_beds = table2array(readtable(append('results/timeseries/hospbeds_scenario', ...
        scenario, '_', tol, 'tol_', date, '.csv')));
    %     TTIQeff_hist = table2array(readtable(append('results/timeseries/TTIQeffhistory_scenario', ...
    %         scenario, '_', tol, 'tol_', date, '.csv')));
    TTIQeff_hist = [];
    AL_hist = table2array(readtable(append('results/timeseries/ALhistory_scenario', ...
        scenario, '_', tol, 'tol_', date, '.csv')));
    
    plotSimLeakyTimeseries2(AL_hist, daily_cases, daily_hosp, cumul_deaths, hosp_beds, TTIQeff_hist, scenario, tol, border)
        
end




