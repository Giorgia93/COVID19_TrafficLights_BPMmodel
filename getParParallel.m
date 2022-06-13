function par = getParParallel(HR_bool, iisol_inf, iborder_inf, ip_trace, itrace_capacity, iVE2, ...
    ireltrans_AL, iptest_clin, iAL0, icomm_seeds, imin_vax_age, icov, ...
    iborder_cases, itol_level, C, tmp, popSizeData)


%------------- Seed and control Parameters --------------
par.R0 = 6; 
par.tEnd = 360;
par.date0 = datenum('01JAN2022');

par.relTransAL = ireltrans_AL; % [0.9, 0.8, 0.7, 0.4] 
par.trafficLights = ["GREEN", "YELLOW", "RED", "EMERGENCY"];
par.AL = iAL0;
par.relTransBaseAL = par.relTransAL(par.AL)*ones(1, par.tEnd+1);

par.minDetectTime = 0;    % time of outbreak detection 
par.followUpTime = 7;     % cases with an isolation time prior to detection are distributed over this time period post-detection
par.borderSeedPeriod = par.tEnd;
par.communitySeedPeriod = 14;
par.borderCases = iborder_cases;
par.meanSeedsPerDay = par.borderCases/365;
par.relInfSeedCases = iborder_inf; %0.2;
par.genSeedVaxStatus = @genSeedVaxStatus_fullyVaxed;

% Active seed cases on day 1
par.nSeeds0 = icomm_seeds;

%%%% Cases tolerance and tier up/down triggers
par.tol_levels = ["low", "medium", "high"];
par.TL = itol_level;

par.trigUp_cases = [50, 200, inf; 100, 400, inf; 300, 1200, inf];
par.trigUp_hospbeds = [inf, inf, 100; inf, inf, 200; inf, inf, 600];
par.trigDown_cases = [0, 100, inf; 75, 300, inf; 200, 800, inf];
par.trigDown_hospbeds = [inf, inf, 50; inf, inf, 150; inf, inf, 400];

% par.trigUp_cases = [10, 200, inf; 25, 400, inf; 50, 1200, inf];
% par.trigUp_hospbeds = [inf, inf, 100; inf, inf, 200; inf, inf, 600];
% par.trigDown_cases = [0, 100, inf; 10, 300, inf; 30, 800, inf];
% par.trigDown_hospbeds = [inf, inf, 50; inf, inf, 150; inf, inf, 400];



%------------- Branching Process Parameters and delays --------------

par.cSub = 0.5; % Relative infectiousness of subclinicals
par.ssk = 0.5; % Overdispersion/superspreading parameter k
par.maxInfectTime = 21; % Maximum infectious period (for computational efficiency)
par.genA = 5.665; par.genB = 2.826; % Generation time distribution parameters
par.incA = 5.8; par.incB = 0.95; % Exposure to Onset distribution parameters
par.isolA = 1; par.isolB = 4; % Onset to isolation distribution parameters
par.traceA = 3; par.traceB = 2/3; % Parent isolation to contact quarantine distribution parameters
par.hospA = 1; par.hospB = 5; % Onset to hospitalisation distribution parameters
par.losA = 1; par.losB = 8; % Hospital LOS distribution parameters

par.pTestClin = iptest_clin; 
par.pTestSub = 0; % Probability of detecting symptomatic / subclinical case
par.pTrace = ip_trace; % Probability of detecting a case by contact tracing
par.traceCapacity = itrace_capacity;

par.cIsol = iisol_inf; %0; %0
par.cQuar = 0.5;


%------------- Vaccine effectiveness --------------
par.VEi1 = 0.55;
par.VEt1 = 0; % Default vaccine efficacy against transmission given infection
par.VEd1 = 0.6; % Default vaccine efficacy against severe disease (and death) given infection
par.VEi2 = iVE2(1); % Default vaccine efficacy against infection
par.VEt2 = iVE2(2); % Default vaccine efficacy against transmission given infection
par.VEd2 = iVE2(3); % Default vaccine efficacy against severe disease (and death) given infection

% Overall effectiveness parameters
par.VEit1 = 1 - (1 - par.VEi1) * (1 - par.VEt1); %0.55
par.VEid1 = 1 - (1 - par.VEi1) * (1 - par.VEd1); %0.82
par.VEit2 = 1 - (1 - par.VEi2) * (1 - par.VEt2); %0.85
par.VEid2 = 1 - (1 - par.VEi2) * (1 - par.VEd2); %0.94


%------------- Disease Rate Data --------------
par.nAgeGroups = 16;
par.IDR = [0.5440, 0.5550, 0.5770, 0.5985, 0.6195, 0.6395, 0.6585, 0.6770, 0.6950, 0.7117, 0.7272, 0.7418, 0.7552, 0.7680, 0.7800, 0.8008]'; % Fraser group
% [par.IHR, par.pICU, par.IFR] = getVerityRates();
[par.IHR, par.pICU, par.IFR] = getHerreraRates();

if HR_bool == true
    OR_tested = 2;
    par.IHR = OR_tested*par.IHR ./ (1-par.IHR+OR_tested*par.IHR);
end

par.ui = [0.4000, 0.3950, 0.3850, 0.4825, 0.6875, 0.8075, 0.8425, 0.8450, 0.8150, 0.8050, 0.8150, 0.8350, 0.8650, 0.8450, 0.7750, 0.7400 ];    % Davies relative susceptibility
par.pSub = 1 - par.IDR;

% fnh = 'pHospFromData_30SEP2021.xlsx';
% fprintf('   Loading hospitalisation rate data:    %s\n', fnh)
% tt = readtable(fnh);
% par.IHR = tt.pHosp;

%------------- Specify Population Structure --------------


par.popCount = zeros(par.nAgeGroups, 1); % Create popDist vector
par.popCount(1:par.nAgeGroups-1) = popSizeData(1:par.nAgeGroups-1, 2); % Fill entries with population distribution
par.popCount(par.nAgeGroups) = sum(popSizeData(par.nAgeGroups:end, 2)); % Aggregate 75+ age-groups
par.totalPopSize = sum(par.popCount);

% get eligibility ratios per age group
if imin_vax_age == 5
    elig_ratios = [0, ones(1, 15)];
elseif imin_vax_age == 12
    elig_ratios = [0, 0, 3/5, ones(1, 13)];
else
    disp("Minimum vax age not valid. Please choose either 5 or 12")
end

par.eligPopCount = par.popCount .* elig_ratios';
% par.eligPopCount = par.popCount .* [0, 0, 3/5, ones(1, 13)]';
%par.totalPopSize = 5.123e6; % Total size of assumed population
par.popDist = par.popCount/sum(par.popCount); 
par.eligPopDist = par.eligPopCount/sum(par.eligPopCount);

par.age = (2.5:5:77.5)'; % Define final age groups (matching contact matrix)


%------------- Define vaccine coverage either static or time-dependent ---------

% Time-dependent scenarios (NOT USED)
% Creates matrices for 1st and 2nd dose coverage by time and age group, starting
% on date0
% If the simulation runs longer than the number of rows in these matrices, the
% coverage on the last defined time date will be used for the remainder of 
% the simulation


% Static scenarios
par.cov1 = zeros(1, par.nAgeGroups);
% par.covElig = 0.90; % proportion of eligible population in the eligible age groups (>=12, 3/5 of 10-14 age group) for 2 doses of vaccine
par.curr_cov2 = [0.0000 0.0000 0.5239 0.8519 0.7981 0.8005 0.8785 0.8966 ...
    0.9212 0.8786 0.9122 0.9028 0.9339 0.9440 0.9605 0.9599]; %[0 0 par.covElig*3/5 par.covElig*ones(1, 13)];

par.cov2 = get_cov(par.curr_cov2, par.eligPopCount, icov, 0);


%------------- Load Contact Matrix and Define NGM --------------
% fs = 'nzcontmatrix.xlsx';
% % fprintf('   Loading contact matrix:    %s\n', fs)
% C = readmatrix(fs); % Get Prem et al contact matrix from data folder
% 
% fs = 'nzpopdist.xlsx';      % This should *ALWAYS* be the national population distribution 'nzpopdist.xlsx'
% % fprintf('   Loading benchmark population distribution for contact matrix:    %s\n', fs)
% tmp = readmatrix(fs); % Load NZ population structure from data folder
popDistBench = [tmp(1:par.nAgeGroups-1, 2); sum( tmp(par.nAgeGroups:end, 2))];

par.C = zeros(par.nAgeGroups, par.nAgeGroups); % Create our contact matrix
for ii = 1:par.nAgeGroups
    for jj = 1:par.nAgeGroups
        par.C(ii,jj) = 0.5*(C(ii,jj) + (popDistBench(jj)/popDistBench(ii)) * C(jj,ii)); % Force detailed balance condition
    end
end

% Re-weight columns for actual population distribution
par.C = (par.popDist./popDistBench).' .* (sum(par.C, 2)./sum((par.popDist./popDistBench).' .* par.C, 2)) .* par.C;

%NGM = par.tI*diag(par.ui)*par.C'*diag(par.IDR + par.cSub*(1-par.IDR)); % Construct "first guess at the NGM"
NGM = diag(par.ui)*par.C'*diag(par.IDR + par.cSub*(1-par.IDR)); % Construct "first guess at the NGM"
par.u = par.R0/max(abs(eig(NGM))); % Choose u to give desired R0
par.NGM = par.u*NGM; % Set final NGM
%par.NGMclin = par.u * par.tI*diag(par.ui)*par.C';   % NGM for clinical individuals (as used in BPM)
par.NGMclin = par.u * diag(par.ui)*par.C';   % NGM for clinical individuals (as used in BPM)


% Proportions of new infections in each group (and by coverage at a given
% time) that single or double vaccinated
pv1 = ((1-par.VEi1)*par.cov1 ./ (1-par.VEi1*par.cov1-par.VEi2*par.cov2));
pv2 = ((1-par.VEi2)*par.cov2 ./ (1-par.VEi1*par.cov1-par.VEi2*par.cov2));

% vaccinated NGM at start date
NGMv0 = (1-par.VEi1*par.cov1(1, :) - par.VEi2*par.cov2(1, :)).' .* (1- par.VEt1*pv1(1, :) - par.VEt2*pv2(1, :)) .* par.NGM;
l = eigs(NGMv0);
Rv = abs(l(1));

end

