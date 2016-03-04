%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization of volcano plot based filering parameters
% -log10(p-value) and log2(fold ratio change) used to analyse random matrix for testing overfitting.
%
% Authors: Peter Horvatovich, Vikram Mitra and Gooitzen Zwanenburg
% Date: February 16, 2016
%
% set mainDir variable to exact location of the ASCA main directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
warning off
% change the directory according to the exact coce location
mainDir = uigetdir;
cd(mainDir)
nRun = 10; % number of repetation

f = importdata('FMatrix.txt'); % experiment design file (the same as the original fractional factorial design)
f([4 11 14],:) = [];
f = f([2,7,5,14,1,15,3,6,12,10,11,8,13,16,4,9],:);
nParam = 2559; % number of variables
sSize = 16; % number of measurements
nFact = 7; % number of factors
const = 6.043; % constant for simulation and for replace 0 in the real data
sd = 0.5391; % standard deviation for simulation
beta = 1; % parameter to calculate f-score

% test ASCA with Fixed pvalue and Fold change thresholds
nRatio = 10;
nPValue = 10;
th = zeros(nRatio,2);
pth = zeros(nPValue);
Dstat = cell(nPValue,nRatio);
sleectedPeaks = zeros(nPValue,nRatio);
simulation = 1; % this parameter is used to optimize parameter for simulated or real experimental data. 0 is real data, 1 is simulated
visualisation = 0; % parameter to allow (1) and disable visualisation (0)
volcanoSelection = 1;% 0 there is no Volcano plot based selection, on the basis of 1 there is a selection
FontSize = 14;% Fontsize in the plots
maxRatioReduction = 40;%reduction of maximal ratio threshold in %
maxPValueReduction = 55;% reduction of maximal p-value threhold in %
mzFilter = 0; %was 400
rtFilterLow = 45; %was 50
rtFilterHigh = 135;%was 130

% setting parallel execution
nPool = 7; % number of threads
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool('local',nPool);
end

parfor kkk = 8:8+7,
    if simulation,
        PvalueTable = cell(1,nRun);
        peaksAffectedNos = zeros(1,nRun);
        R_main = zeros(sSize ,nParam);
        peaksAffected = zeros(nFact,nParam);
        factorLabel = {'Factor 1','Factor 2','Factor 3','Factor 4','Factor 5','Factor 6 ','Factor 7'};
        
        for mn = 1:nFact,
            R_main = R_main+normrnd(0,sd,sSize,nParam);
        end
        
        R_main = R_main/sqrt(nFact);
        R_main = R_main + const;
        R_main_10p=exp(R_main);
        
        if find(R_main == 0),
            disp('0 in the matrix');
        end;
    else
        % importing data from the mpks file
        load('factorLabel.mat')
        factorLabel{2} = ['\bf' factorLabel{2}];
        factorLabel{5} = ['\it' factorLabel{5}];
        factorLabel{6} = ['\bf' factorLabel{6}];
        %     for i = 1:length(factorLabel),
        %         factorLabel{i} = [factorLabel{i} ' '];
        %     end
        data = importdata('Outstem_mzRadius=0.3_TRadius=1_Fraction=0.5_mzStart=100_rtStart=65_mzEnd=1500_rtEnd=135.mpks');
        %     data = importdata('d:\users\peter\munka\publications\2012\vikram_fact_des\code\codeOk\Outstem_mzRadius=0.3_TRadius=1_Fraction=0.5_mzStart=100_rtStart=65_mzEnd=1500_rtEnd=135_16f_ok.mpks');
        % removing data points < 400da Mz and < 50Mins and >130Mins
        data_FS = data.data;
        
        %	data_FS(:,[4 11 14])=[]; % use only with data including 4 repepatition
        [filtdata_mz indices] = find(data_FS(:,2) > mzFilter);
        filtdata_mz = data_FS(filtdata_mz,:);
        [filtdata_rt indices] = find(filtdata_mz(:,3) > rtFilterLow & filtdata_mz(:,3) <rtFilterHigh);
        filtdata_rt = filtdata_mz(filtdata_rt,:);
        dat_FS = filtdata_rt(:,14:end)';
        % resampling and replacing missing intensity values
        NonZeroMat = dat_FS;
        for i = 1:size(NonZeroMat,2)
            zeroIn = find(NonZeroMat(:,i) == 0);
            nZeros = length(zeroIn);
            R = lognrnd(const,sd,[size(NonZeroMat,1),1]);
            NonZeroMat(zeroIn,i) = R(zeroIn);
        end
        R_main_10p=NonZeroMat;
        peaksAffected = [];
    end
    
    % setup P-value and log 2 ratio threholds
    [maxPRatio] = volcano_plot(R_main_10p,f,0,0,0,[],1,0,factorLabel);
    th1 = linspace(maxPRatio(2)*(100-maxRatioReduction)/100,0,nRatio)'; % log2 fold + ratio threhold
    th2 = -th1; % log2 fold - ratio threhold
    th = [th2,th1]; % log2 fold ratio threhold
    pth = linspace(maxPRatio(1)*(100-maxPValueReduction)/100,0,nPValue)'; % log10 p-value threhold
    
    loopParamOpt(th, pth, R_main_10p, f, nFact, visualisation, factorLabel, nParam, kkk, sSize, simulation, FontSize, nRatio, nPValue, mainDir);
end
delete(poolobj)
