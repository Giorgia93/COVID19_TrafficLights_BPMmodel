function [nInfected, nCases, nHosp, nHospBeds, nDeaths, ...
    AL_timeseries, AL_props, TTIQeff, TTIQeff_time] = runSimLeaky_LessVerbose(par, earlyReject, verbose)


% transmission process will step if this many people have been infected
maxCases = 1000000;

% Set up time array
t = par.date0 + (0:1:par.tEnd);

% Set the relative transmission rate due to alert level settings, as a function of time
% (this is the base setting and can be later modified by dynamic control rules)
relTransCurrentAL = par.relTransBaseAL;

% calculate the area under the curve of the generation time distribution in 1-day time
% steps - this gives the relative amount of transmitting each person does
% on each day of their infection
tArr = 0:1:par.maxInfectTime+1;
C = wblcdf(tArr, par.genA, par.genB);
auc = diff(C).';

% initialise variables for the case table
caseID = (1:maxCases).';
parentID = nan(maxCases, 1);
gen = nan(maxCases, 1);
nOff = zeros(maxCases, 1);
ageGroup = nan(maxCases, 1);
Rimult = genRi(maxCases, 1, par);   % pre generate each individual's values of Ri
subclinFlag = nan(maxCases, 1);
vaccDoses = zeros(maxCases, 1);
tInfect = nan(maxCases, 1);
tOnset = genOnsetDelay(maxCases, 1, par);   % pre-generate each individuals incubation period
tIsol = nan(maxCases, 1);   % testing and isolation times assumed simulteneous
tQuar = nan(maxCases, 1);   % testing and isolation times assumed simulteneous
tHosp = nan(maxCases, 1);
tDisc = nan(maxCases, 1); 
icuFlag = zeros(maxCases, 1); 
diedFlag = zeros(maxCases, 1); 

% initialise case table 
cases = table(caseID, parentID, gen, nOff, ageGroup, Rimult, subclinFlag, vaccDoses, tInfect, tOnset, tQuar, tIsol, tHosp, tDisc, icuFlag, diedFlag);

%%%% define number of seed cases and their properties

% Draw total number of seed cases over the whole simulated period
nBorderCases = poissrnd(par.meanSeedsPerDay*par.borderSeedPeriod) + 0;

% Total number of seed cases is the existing seed cases at day 1 + border
% cases arriving during the year
nSeedCases = par.nSeeds0 + nBorderCases;

% Set parentID of all seed cases to 0
cases.parentID(1:nSeedCases) = 0;

% Set generation of seed cases to 1 
cases.gen(1:nSeedCases) = 1;

% Set ages of seed cases to be distributed the same as the overall popn
cases.ageGroup(1:nSeedCases) = 1 + discRand(par.popDist.', nSeedCases, 1);    

% Reduce infectiousness of border cases by multiplying Rimult by 
% relInfSeedCases <1, to represent home isolation of border cases
cases.Rimult(par.nSeeds0+1:nSeedCases) = par.relInfSeedCases * cases.Rimult(par.nSeeds0+1:nSeedCases); 

% Set subclinical status for seed cases from probability of being
% subclinical for each age group
cases.subclinFlag(1:nSeedCases) = rand(nSeedCases, 1) < par.pSub(cases.ageGroup(1:nSeedCases));
% Set vaccination status of seed cases
seed0_cov = par.cov2 .* (1 - par.VEi2) ./ (1 - (par.cov2 .* par.VEi2));
cases.vaccDoses(1:par.nSeeds0) = 2 .* (rand(1, par.nSeeds0) < seed0_cov(cases.ageGroup(1:par.nSeeds0)));
cases.vaccDoses(par.nSeeds0+1:nSeedCases) = par.genSeedVaxStatus(cases, nBorderCases, par);

% Set time of appearance of each seed case, randomly drawn from the
% simulated time period
cases.tInfect(1:par.nSeeds0) = t(1) + floor(par.communitySeedPeriod*rand(par.nSeeds0, 1));
cases.tInfect(par.nSeeds0+1:nSeedCases) = t(1) + floor(par.borderSeedPeriod*rand(nBorderCases, 1));

% Set incubation/generation (?) period of seed cases from randomly drawn inc periods
cases.tOnset(1:nSeedCases) = cases.tOnset(1:nSeedCases) + cases.tInfect(1:nSeedCases);

% Set probability that seed cases go into isolation from the probabilities 
% of symptomatic and subclinical cases getting tested
pIsol = par.pTestClin*(cases.subclinFlag(1:nSeedCases) == 0 ) + par.pTestSub*(cases.subclinFlag(1:nSeedCases) == 1 );

% Set isolation status of seed cases 
isolFlag = rand(nSeedCases, 1) < pIsol;
nIsol = sum(isolFlag);

% Set time when seed cases went into isolation
tIsol = cases.tOnset(isolFlag) + genIsolDelay(nIsol, 1, par);

% Find which seed cases started isolating before the minimum detection time
ind = tIsol < t(1)+par.minDetectTime;
% Make sure all seed cases start isolating after min detection time
tIsol(ind) = t(1)+par.minDetectTime + rand(sum(ind), 1)*par.followUpTime;
cases.tIsol(isolFlag) = tIsol;

% Fractions of each age group that is in the susceptible class and has had
% 0, 1 or 2 doses of vaccines. These will initially sum to 1 but will
% become depleted as people get infected
susFrac0 = 1-par.cov1(1, :)-par.cov2(1, :);
susFrac1 = par.cov1(1, :);
susFrac2 = par.cov2(1, :);

% Fractions moving between dose compartments each time step
q01 = min(1, max(0, (par.cov1(2:end, :)+par.cov2(2:end, :)-par.cov1(1:end-1, :)-par.cov2(1:end-1, :))./(1-par.cov1(1:end-1, :)-par.cov2(1:end-1, :))));
q01 = [q01; zeros(1, par.nAgeGroups)];
q12 = min(1, max(0, (par.cov2(2:end, :)-par.cov2(1:end-1, :))./par.cov1(1:end-1, :)));
q12 = [q12; zeros(1, par.nAgeGroups)];
[nRows_qij, ~] = size(q01);

% Relative transmisison rates as a result of vaccination and
% quarantine/isolation
relTransVacc = [1; 1-par.VEt1; 1-par.VEt2];
relTransIsol = [1; par.cQuar; par.cIsol];


nSteps = length(t)-1;

% Number of seed cases in each age group
nCases0 = histcounts(cases.ageGroup(1:nSeedCases), 1:par.nAgeGroups+1);
nCases = nCases0;
nActive = 0;
% nActiveDetected = zeros(size(t));
nDailyAvg = zeros(size(t));
nHospBeds = zeros(size(t));
TTIQeff_time = zeros(size(t));
nFuture = nSeedCases;
iStep = 1;
dist = 0;


% Simulation stops when:
% - Max number of steps reached OR
% - Max number of cases reached OR
% - Total number of seed cases over entire simulated time = 0 and 
%   the relative transmission at the current AL = 1 OR
% - The sum of squares btw our sim results and data goes beyond set threshold

while iStep < nSteps && sum(nCases) < maxCases && (nActive + nFuture > 0 || relTransCurrentAL(iStep+1) < 1) && dist < earlyReject.threshold 
    
    
    iStep = iStep+1;
    
    % number of seed cases whose infeciotn time is still in the future:
    nFuture = sum( t(iStep) <= cases.tInfect );
    
    % Get IDs of active cases at current time step:
    activeID = cases.caseID(t(iStep) > cases.tInfect & t(iStep) <= cases.tInfect+par.maxInfectTime);
    nActive = length(activeID);
    
    if nActive > 0
        % Simulate new secondary cases on day i from each currently active
        % case:
        
        % Area under the curve of the transmission rate for each active case at current time step
        auci = auc(t(iStep)-cases.tInfect(activeID));       
        
        % isolation status of active cases (0 for nothing, 1 for
        % quarantined, 2 for isolated)
        isolStatus = (t(iStep) > cases.tQuar(activeID) & ~(t(iStep) > cases.tIsol(activeID))) + 2*(t(iStep) > cases.tIsol(activeID));
        
        % Calculate the expected number of offspring in each age group during the current
        % time step from each active case, assuming a fully susceptible population.
        % This defines a matrix expOff whose i,j element is the expected number of
        % offspring from parent case i in age group j
        % This is the product of the following factors for each case:
        % - relative transmission due to current alert level
        % - relative transmission due to vaccination status of active case i
        % - relative reproduction number (typically gamma distributed) of case i
        % - relative transmission rate for clinical/subclinical status of case i
        % - relative transmission due to isolation/quarantine status of case i
        % - time-dependent transmission rate for the parent case i, quantified by auci
        % - jth column of the NGM for an unvaccinated clinical individual in the age group of the parent case
        expOff = relTransCurrentAL(iStep) * relTransVacc(1+cases.vaccDoses(activeID)) .* relTransIsol(1+isolStatus) .* cases.Rimult(activeID) .* (1 - (1-par.cSub)*cases.subclinFlag(activeID)) .* auci .* (par.NGMclin(:, cases.ageGroup(activeID))).';
        
        % split expOff into expected offspring who have had 0, 1 or 2 doses
        % In doing this, some putative infections are prevented by the
        % vaccine (and others by immunity from prior infection is susFrac0,
        % susFrac1 and susFrac2 sum to less than 1).
        expOff0 = expOff .* susFrac0;
        expOff1 = (1-par.VEi1) * expOff .* susFrac1;
        expOff2 = (1-par.VEi2) * expOff .* susFrac2;
        
        % Generate nOff actual number of offspring from parent case i in age group j, 
        % thinned due to immunity in each age group as a result of infection prevention 
        % in vaccinated individuals and prior infection prevention in vaccinated offspring
        nOff0 = poissrnd(expOff0);   
        nOff1 = poissrnd(expOff1); 
        nOff2 = poissrnd(expOff2); 
        
        % Total number of offspring summed across all parent cases and all
        % age groups:
        nOffTot = sum(sum(nOff0+nOff1+nOff2));
        
        % If number of new cases exceeds maxCases, reduce them to prevent table overflow 
        if nOffTot > maxCases-sum(nCases)     
            excess = nOffTot - (maxCases-sum(nCases));
            M = [nOff0(:); nOff1(:); nOff2(:)];
            X = repelem((1:length(M)).', M, 1);
            deletionIndices = X(randsample(length(X), excess));
            deletions = histcounts(deletionIndices, 1:length(M)+1).';
            M = M-deletions;
            nOff0 = reshape(M(1:(nActive*par.nAgeGroups)), nActive, par.nAgeGroups);
            nOff1 = reshape(M((nActive*par.nAgeGroups+1):(2*nActive*par.nAgeGroups)), nActive, par.nAgeGroups);
            nOff2 = reshape(M((2*nActive*par.nAgeGroups+1):end), nActive, par.nAgeGroups);
            nOffTot = sum(sum(nOff0+nOff1+nOff2));
        end
        
        if nOffTot > 0
            secIDs = (sum(nCases)+1:sum(nCases)+nOffTot).';       % IDs for today's newly infected cases

            % assign age groups, parent ID, clinical status and vaccination status for new
            % cases based on nOff matrices
            ageGroupList = repmat(1:par.nAgeGroups, nActive, 1);
            parentNumber = repmat(1:nActive, 1, par.nAgeGroups);
            nx = par.nAgeGroups*nActive;
            % Make an array (propList) whose coluumns are: (1) age group of
            % offspring, (2) ID of parent case, (3) number of doses of offspring
            propList = [repmat([ageGroupList(:), parentNumber(:)], 3, 1), [zeros(nx, 1); ones(nx, 1); 2*ones(nx, 1)] ];
            % The frequency (number of cases with each age gorup, parent ID
            % and dose combination) of each row of propList is given by
            % the elements in the nOff matrices. Generate a sample X whose
            % columns are the values of these three properties for each
            % secondary case:
            X = repelem( propList, [nOff0(:); nOff1(:); nOff2(:)], 1);
            cases.ageGroup(secIDs) = X(:, 1);
            cases.parentID(secIDs) = activeID( X(:, 2) );
            cases.vaccDoses(secIDs) = X(:, 3);
            cases.subclinFlag(secIDs) = rand(nOffTot, 1) < par.pSub(cases.ageGroup(secIDs));
            
            cases.gen(secIDs) = cases.gen(cases.parentID(secIDs))+1;    % Generation of each new cases is generation of parent + 1
            cases.tInfect(secIDs) = t(iStep);                           % Infection time for each new cases is today
            cases.tOnset(secIDs) = t(iStep) + cases.tOnset(secIDs);     % For efficiency infection to onset delay is pre-stored in cases.tOnset
            cases.nOff(activeID) = cases.nOff(activeID)+sum(nOff0+nOff1+nOff2, 2);      

            % simulate case testing and isolation effects for new cases
            pIsol = par.pTestClin*(cases.subclinFlag(secIDs) == 0 ) + par.pTestSub*(cases.subclinFlag(secIDs) == 1 ) ;
            isolFlag = rand(nOffTot, 1) < pIsol;
            nIsol = sum(isolFlag);
            tIsol = cases.tOnset(secIDs(isolFlag)) + genIsolDelay(nIsol, 1, par);
            ind = tIsol < t(1)+par.minDetectTime;
            tIsol(ind) = t(1)+par.minDetectTime + rand(sum(ind), 1)*par.followUpTime;
            cases.tIsol(secIDs(isolFlag)) = tIsol;
            
            % simulate contact tracing effects for new cases
%             pTrace = par.pTrace * (~(nActiveDetected(iStep-1) > par.traceCapacity)) *(~isnan(cases.tIsol(cases.parentID(secIDs))));
            pTrace = par.pTrace * (~(nDailyAvg(iStep-1) > par.traceCapacity)) *(~isnan(cases.tIsol(cases.parentID(secIDs))));
            traceFlag = rand(nOffTot, 1) < pTrace;
            nTrace = sum(traceFlag);
            cases.tQuar(secIDs(traceFlag)) = cases.tIsol(cases.parentID(secIDs(traceFlag))) + genTraceDelay(nTrace, 1, par);
            % optional: individuals traced prior to onset go into full isolation (as opposed to quarantine) on symptom onset:
            % this is applied to all individuals including subclinical on
            % assumption that asymtomatic contacts get tested. This allows
            % offspring of subclinicals to be traced
            cases.tIsol(secIDs(traceFlag)) = min(cases.tIsol(secIDs(traceFlag)), max(cases.tQuar(secIDs(traceFlag)), cases.tOnset(secIDs(traceFlag))));

            % simulate clinical outcomes (hsopitalisation, ICU and death) for new cases
            clinFlag = cases.subclinFlag(secIDs) == 0;
            pHosp = ((cases.vaccDoses(secIDs) == 0) + ...
                (1-par.VEd1)*(cases.vaccDoses(secIDs) == 1) + ...
                (1-par.VEd2)*(cases.vaccDoses(secIDs) == 2)) .* par.IHR(cases.ageGroup(secIDs))./par.IDR(cases.ageGroup(secIDs));
            hospFlag = clinFlag & (rand(nOffTot, 1) < pHosp);
            cases.icuFlag(secIDs) = hospFlag & (rand(nOffTot, 1) < par.pICU(cases.ageGroup(secIDs))  );  
            cases.diedFlag(secIDs) = hospFlag & (rand(nOffTot, 1) < par.IFR(cases.ageGroup(secIDs))./par.IHR(cases.ageGroup(secIDs)));  
            cases.tHosp(secIDs(hospFlag)) = cases.tOnset(secIDs(hospFlag)) + genHospDelay(sum(hospFlag), 1, par);
            cases.tDisc(secIDs(hospFlag)) = cases.tHosp(secIDs(hospFlag)) + genHospLOS(sum(hospFlag), 1, par);
            % cases are detected and isolated once hospitalised
            cases.tIsol(secIDs(hospFlag)) = min(cases.tIsol(secIDs(hospFlag)), cases.tHosp(secIDs(hospFlag)));
           
            % Update cumulative infections to date (in each age group)
            nCases = nCases+sum(nOff0+nOff1+nOff2, 1);        
        end
    else
        nOff0 = zeros(1, par.nAgeGroups);
        nOff1 = zeros(1, par.nAgeGroups);
        nOff2 = zeros(1, par.nAgeGroups);
    end   

    % Update susceptible fractions according to vaccinations given
    % and new infections this time step
    iRow = min(iStep-1, nRows_qij);
    susFrac2 = susFrac2 + susFrac1.*q12(iRow, :) - sum(nOff2, 1)./(par.popCount');
    susFrac1 = susFrac1.*(1-q12(iRow, :)) + susFrac0.*q01(iRow, :) - sum(nOff1, 1)./(par.popCount');
    susFrac0 = susFrac0.*(1-q01(iRow, :)) - sum(nOff0, 1)./(par.popCount');
   
    
    
    %%%%%% Check triggers and update current alert level accordingly %%%%%%

    % Average number of new daily cases over the last 7 days
    nDailyAvg(iStep) = sum( t(iStep) > min(cases.tIsol, cases.tQuar) & t(iStep) <= min(cases.tIsol, cases.tQuar)+7 ) / 7;
    % Current number of occupied hospital beds
    nHospBeds(iStep) = sum(t(iStep) > cases.tHosp & t(iStep) <= cases.tDisc);
    
    % Move up levels if trigger met (stay at initial AL for first 28 days)
    if par.AL < 4 && iStep > 28 && (nDailyAvg(iStep) > par.trigUp_cases(par.AL, par.TL) || nHospBeds(iStep) > par.trigUp_hospbeds(par.AL, par.TL))
        par.AL = par.AL + 1;
        relTransCurrentAL(iStep+1:end) = par.relTransAL(par.AL);
        
    % Move down levels if trigger met
    elseif par.AL > 1 && iStep > 28 && (nDailyAvg(iStep) <= par.trigDown_cases(par.AL - 1, par.TL) > 0 && nHospBeds(iStep) <= par.trigDown_hospbeds(par.AL - 1, par.TL))
        par.AL = par.AL - 1;
        relTransCurrentAL(iStep+1:end) = par.relTransAL(par.AL);
    end
    
    
    
    % calculate distance metric to allow early rejection of simulation
    % during fitting process - NOT CURRENTLY USED
    nIsolTemp = histcounts(cases.tIsol, [t(1):t(iStep)+1] );
    dist = calcError(earlyReject.tData, earlyReject.nCasesData, t(1):t(iStep), nIsolTemp);
    


end
nActiveEnd = nActive;

if iStep < nSteps
   dist = 2*earlyReject.threshold; 
end

totCases = sum(nCases);
cases = cases(1:totCases, :);

% Only count Reff over cases who have completed their infectious period
Reff = mean(cases.nOff(cases.tInfect < t(iStep) - par.maxInfectTime));

% % calculate TTIQ effect overall and over time
% wI = 1-wblcdf(floor(cases.tIsol)-cases.tInfect, par.genA, par.genB); % amount of transmission prevented if isolation was perfect
% wQ = 1-wblcdf(floor(min(cases.tQuar, cases.tIsol))-cases.tInfect, par.genA, par.genB); % amount of transmission prevented if quarantine was perfect    % wQ >= wI
% wI(isnan(wI)) = 0;
% wQ(isnan(wQ)) = 0;
% ReffReduc = (1-par.cQuar)*(wQ-wI) + (1-par.cIsol)*(wI);
% infAfterDetectFlag = cases.tInfect >= t(1)+par.minDetectTime;
% pClin = mean(cases.subclinFlag(infAfterDetectFlag) == 0);
% ReffReducClin = mean(ReffReduc(infAfterDetectFlag & cases.subclinFlag == 0));
% ReffReducSub = mean(ReffReduc(infAfterDetectFlag & cases.subclinFlag == 1));
% TTIQeff = 1 - (pClin*(1-ReffReducClin) + (1-pClin)*par.cSub*(1-ReffReducSub)) / (pClin + par.cSub*(1-pClin));
% infUnderCapFlag = interp1(t, nDailyAvg, max(0, cases.tInfect-1)) <= par.traceCapacity;
% ReffReducClin = mean(ReffReduc(infAfterDetectFlag & infUnderCapFlag & cases.subclinFlag == 0));
% ReffReducSub = mean(ReffReduc(infAfterDetectFlag & infUnderCapFlag & cases.subclinFlag == 1));
% TTIQeffUnderCap = 1 - (pClin*(1-ReffReducClin) + (1-pClin)*par.cSub*(1-ReffReducSub)) / (pClin + par.cSub*(1-pClin));
% ReffReducClin = mean(ReffReduc(infAfterDetectFlag & ~infUnderCapFlag & cases.subclinFlag == 0));
% ReffReducSub = mean(ReffReduc(infAfterDetectFlag & ~infUnderCapFlag & cases.subclinFlag == 1));
% TTIQeffOverCap = 1 - (pClin*(1-ReffReducClin) + (1-pClin)*par.cSub*(1-ReffReducSub)) / (pClin + par.cSub*(1-pClin));
% TTIQeff_time(nDailyAvg <= par.traceCapacity) =  TTIQeffUnderCap;
% TTIQeff_time(nDailyAvg > par.traceCapacity) = TTIQeffOverCap;

% calculate TTIQ effect overall and over time
wI = 1-wblcdf(floor(cases.tIsol)-cases.tInfect, par.genA, par.genB); % amount of transmission prevented if isolation was perfect
wQ = 1-wblcdf(floor(min(cases.tQuar, cases.tIsol))-cases.tInfect, par.genA, par.genB); % amount of transmission prevented if quarantine was perfect    % wQ >= wI
wI(isnan(wI)) = 0;
wQ(isnan(wQ)) = 0;
ReffReduc = (1-par.cQuar)*(wQ-wI) + (1-par.cIsol)*(wI);
infAfterDetectFlag = cases.tInfect >= t(1)+par.minDetectTime;
indSub = infAfterDetectFlag & cases.subclinFlag == 0;
indClin = infAfterDetectFlag & cases.subclinFlag == 1;
TTIQeff = 1 - (sum(1-ReffReduc(indClin)) + par.cSub*sum(1-ReffReduc(indSub))) / (sum(indClin) + par.cSub*sum(indSub));
TTIQeff_time = zeros(size(t));

for ii = 1:length(t)
    indSub = infAfterDetectFlag & cases.subclinFlag == 0 & (cases.tInfect <= t(ii) & cases.tInfect > t(ii)-7);
    indClin = infAfterDetectFlag & cases.subclinFlag == 1 & (cases.tInfect <= t(ii) & cases.tInfect > t(ii)-7);
    TTIQeff_time(ii) = 1 - (sum(1-ReffReduc(indClin)) + par.cSub*sum(1-ReffReduc(indSub))) / (sum(indClin) + par.cSub*sum(indSub));
end

TTIQeff_time(isnan(TTIQeff_time)) = 0;


% Count number of cases (and other outcomes) in each age group on each day
tExt = [t t(end) + 1];
nInfected = histcounts2(cases.tInfect, cases.ageGroup, tExt, 1:par.nAgeGroups+1);
nIsol = histcounts2(cases.tIsol, cases.ageGroup, tExt, 1:par.nAgeGroups+1);
nHosp = histcounts2(cases.tHosp, cases.ageGroup, tExt, 1:par.nAgeGroups+1);
nDisc = histcounts2(cases.tDisc, cases.ageGroup, tExt, 1:par.nAgeGroups+1);
nICUIn = histcounts2(cases.tHosp(cases.icuFlag == 1), cases.ageGroup(cases.icuFlag == 1), tExt, 1:par.nAgeGroups+1);
nICUOut = histcounts2(cases.tDisc(cases.icuFlag == 1), cases.ageGroup(cases.icuFlag == 1), tExt, 1:par.nAgeGroups+1);
nDeaths = histcounts2(cases.tDisc(cases.diedFlag == 1), cases.ageGroup(cases.diedFlag == 1), tExt, 1:par.nAgeGroups+1);

% Daily number of detected cases calculated for each age group as infected
% individual who had a non-null quarantine time or isolation time
nCases = histcounts2(cases.tInfect(~isnan(cases.tQuar) | ~isnan(cases.tIsol)), cases.ageGroup(~isnan(cases.tQuar) | ~isnan(cases.tIsol)), tExt, 1:par.nAgeGroups+1);

% Change in AL over time
AL_timeseries = interp1(par.relTransAL,1:numel(par.relTransAL),relTransCurrentAL);
AL_props = [sum(AL_timeseries == 1)/length(AL_timeseries) ...
    sum(AL_timeseries == 2)/length(AL_timeseries)...
    sum(AL_timeseries == 3)/length(AL_timeseries)...
    sum(AL_timeseries == 4)/length(AL_timeseries)];

if verbose == true
    fprintf("infections cases hosp peakBeds deaths G Y R E TTIQeff\n")
    fprintf("  %i   %i %i  %i   %i    %i%% %i%% %i%% %i%%  %i%%\n", ...
        sum(nInfected, 'all'), sum(nCases, 'all'), ...
        sum(nHosp, 'all'), max(nHospBeds), sum(nDeaths, 'all'), ...
        floor(100*AL_props(1)), floor(100*AL_props(2)), ...
        floor(100*AL_props(3)), floor(100*AL_props(4)), floor(100*TTIQeff))
end

end



