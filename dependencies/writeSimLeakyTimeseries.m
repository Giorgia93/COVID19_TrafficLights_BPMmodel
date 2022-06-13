function [] = writeSimLeakyTimeseries(daily_inf, daily_cases, daily_hosp, cdeaths, hosp_beds, AL_history, TTIQeff_history, scenario, tol, today)

filename_dailyinf = append('results/timeseries/dailyinf_scenario', scenario, '_', ...
    tol, 'tol_', today, '.csv');
writetable(array2table(daily_inf), filename_dailyinf)

filename_dailycases = append('results/timeseries/dailycases_scenario', scenario, '_', ...
    tol, 'tol_', today, '.csv');
writetable(array2table(daily_cases), filename_dailycases)

filename_dailyhosp = append('results/timeseries/dailyhosp_scenario', scenario, '_', ...
    tol, 'tol_', today, '.csv');
writetable(array2table(daily_hosp), filename_dailyhosp)

filename_cdeaths = append('results/timeseries/cumuldeaths_scenario', scenario, '_', ...
    tol, 'tol_', today, '.csv');
writetable(array2table(cdeaths), filename_cdeaths)

filename_hospbeds = append('results/timeseries/hospbeds_scenario', scenario, '_', ...
    tol, 'tol_', today, '.csv');
writetable(array2table(hosp_beds), filename_hospbeds)

filename_ALhistory = append('results/timeseries/ALhistory_scenario', scenario, '_', ...
    tol, 'tol_', today, '.csv');
writetable(array2table(AL_history), filename_ALhistory)

filename_TTIQeff_history = append('results/timeseries/TTIQeffhistory_scenario', scenario, '_', ...
    tol, 'tol_', today, '.csv');
writetable(array2table(TTIQeff_history), filename_TTIQeff_history)

end