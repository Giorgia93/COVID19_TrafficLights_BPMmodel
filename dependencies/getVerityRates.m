function [IHR, pICU, IFR] = getVerityRates()

IHR = [1.0200E-04, 1.0200E-04, 3.0600E-04, 2.9060E-03, 7.9020E-03, 1.6375E-02, 2.8325E-02, 3.6350E-02, 4.0450E-02, 5.2275E-02, 7.1825E-02, 9.0700E-02, 1.0890E-01, 1.3000E-01, 1.5400E-01, 1.7809E-01]';   % Verity (note 1st element need to be nonzero to avoid divide by zero error) - set 0-4 rates equal to 5-9 rates
IFR = [2.9450E-05, 2.9450E-05, 5.6150E-05, 1.2938E-04, 2.4913E-04, 4.4275E-04, 7.1025E-04, 1.0430E-03, 1.4410E-03, 2.7175E-03, 4.8725E-03, 9.2875E-03, 1.5963E-02, 2.5175E-02, 3.6925E-02, 6.6433E-02]';    % Verity 
pICU = [0.11 0.14 0.16 0.18 0.21 0.24 0.28 0.32 0.37 0.41 0.45 0.48 0.45 0.39 0.31 0.20].';     % Probability of ICU given hospitalisation Knock et al

% Apply Delta hazard ratios from Twohig et al & Fisman et al
HR_Hosp_Delta = 2.26;
OR_Death_Delta = 2.32;
IHR = 1 - (1-IHR).^HR_Hosp_Delta;
IFR = OR_Death_Delta * IFR ./ (1 - IFR + OR_Death_Delta * IFR);

end