%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization of volcano plot based filering parameters
% -log10(p-value) and log2(fold ratio change).
% script can be used for simulated and real data (use simulation parameter to
% choose between these two options (0 is real data, 1 is simulated).
% Figure S2-4 was made with this script
%
% Authors: Peter Horvatovich, Vikram Mitra and Gooitzen Zwanenburg
% Date: July 17, 2015
%
% set mainDir variable to exact location of the ASCA main directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
warning off
% change the directory according to the exact coce location
%mainDir = 'e:\data\publications\2012\vikram_fact_des\code\codeOkPerm\';%uigetdir;
% change the directory according to the exact code location
mainDir = uigetdir;
cd(mainDir)
nRun = 1; % number of repetation

f = importdata('FMatrix.txt'); % experiment design file (the same as the original fractional factorial design)
f([4 11 14],:) = [];
f = f([2,7,5,14,1,15,3,6,12,10,11,8,13,16,4,9],:);
nParam = 1000; % number of variables
sSize = 16; % number of measurements
nFact = 7; % number of factors
const = 6.043; % constant for simulation and for replace 0 in the real data
sd = 0.5391; % standard deviation for simulation
beta = 1; % parameter to calculate f-score
lowM = [-1.5,0,2,0,-3,0,0]; % mean for lower level
highM = [1.5,0,-2,0,3,0,0]; % mean for higher level
percentagePeaks = [0.05, 0, 0.035, 0, 0.05, 0, 0]; % % of factor affected peaks
peaksNAffected = round(repmat(nParam, nFact,1).*percentagePeaks')'; % vector contaning the number of factor influenced peaks
peaksAffected = cell(nRun,1); % matrix containing the index of peaks, which are influenced by the factor (with 1)

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

if simulation,
    Recall = zeros(length(pth),length(th));
    precision = zeros(length(pth),length(th));
    contingencyTable = cell(length(pth),length(th));
    specificity = zeros(length(pth),length(th));
    gScore = zeros(length(pth),length(th));
    fScore = zeros(length(pth),length(th));
    PvalueTable = cell(1,nRun);
    peaksAffectedNos = zeros(1,nRun);
    R_main = zeros(sSize ,nParam);
    peaksAffected = zeros(nFact,nParam);
    factorLabel = {'\bfFactor 1','Factor 2','\bfFactor 3','Factor 4','\bfFactor 5','Factor 6 ','Factor 7'};
    
    for nn=1:nFact,
        Dummy=randperm(nParam);
        peaksAffected(nn, Dummy(1,1:peaksNAffected(nn)))=ones(1,peaksNAffected(nn));
    end
    peaksUniqueAffected=max(peaksAffected);
    tindx = find(peaksUniqueAffected);
    nindx = find(peaksUniqueAffected==0);
    
    for mn = 1:nFact,
        % peaks affceted
        high = find(f(:,mn) == 2);
        low = find(f(:,mn) == 1);
        if peaksNAffected(mn)>0,
            R_main(low,find(peaksAffected(mn,:)))=R_main(low,find(peaksAffected(mn,:)))+normrnd(lowM(mn),sd,[size(low,1),peaksNAffected(mn)]);
            R_main(high,find(peaksAffected(mn,:)))= R_main(high,find(peaksAffected(mn,:)))+normrnd(highM(mn),sd,[size(high,1),peaksNAffected(mn)]);
        end
        R_main(:,find(peaksAffected(mn,:)==0))= R_main(:,find(peaksAffected(mn,:)==0))+normrnd(0,sd,sSize,[size(find(peaksAffected(mn,:)==0),2)]);
    end
    
    peaksAffectedNos = length(tindx);
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
[maxPRatio] = volcano_plot(R_main_10p,f,0,0,0,peaksAffected,1,0,factorLabel);
th1 = linspace(maxPRatio(2)*(100-maxRatioReduction)/100,0,nRatio)'; % log2 fold + ratio threhold
th2 = -th1; % log2 fold - ratio threhold
th = [th2,th1]; % log2 fold ratio threhold
pth = linspace(maxPRatio(1)*(100-maxPValueReduction)/100,0,nPValue)'; % log10 p-value threhold

% Test ASCA with variying fold change
for i = 1:length(th),
    for j = 1:length(pth),
        %ASCA with permutation p-value
        [XASCA,indx,p_dat] = ASCAwithPerm(R_main_10p,f,sSize, [1:nFact], pth(j),th(i,:),[],visualisation, factorLabel);
        sleectedPeaks(i,j)=size(indx,1);
        indxNotSel=1:nParam;
        indxNotSel(indx)=[];
        ASCA{i,j} = XASCA;
        D = XASCA.effects.ssq_factors';
        disp(['Iteration: ' num2str(i) ' and ' num2str(j)]);
        if ~(isempty(XASCA)),
            Dstat{i,j} = D(:,1);
            PvalueTable{i,j} = p_dat;
            if simulation,
                % tindx true indeces and nindx are false indeces
                tp = length(intersect(tindx,indx));
                fp = length(intersect(nindx,indx));
                tn = length(intersect(nindx,indxNotSel));
                fn = length(intersect(tindx,indxNotSel));
                contingencyTable{i,j}=[tp,fp;fn,tn];
                if (tp+fp+tn+fn)~=nParam,
                    disp('Problem with contigency table!')
                end
                
                % Calculate specificity and sensitivity
                recall(i,j) = tp/(tp+fn);
                precision(i,j) = tp/(tp+fp);
                specificity(i,j) = tn/(tn+fp);
                gScore(i,j)=sqrt(specificity(i,j)*recall(i,j));
                fScore(i,j)=((1+beta^2)*(precision(i,j)*recall(i,j)))/((beta^2)*precision(i,j)+recall(i,j));
            end
        else
            Dstat{i,j} = zeros(nFact,1);
            recall(i,j) = 0;
            precision(i,j) = 0;
            PvalueTable{i,j} = ones(nFact,1);
        end
        XX(i,j)=th(i,1);
        XX2(i,j)=th(i,2);
        YY(i,j)=pth(j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualisation
% 
% This code part is used to load saved simulated data and make visualisation
% clc
% clear
% load('e:\data\publications\2012\vikram_fact_des\code\codeOkPerm\ParamOptSimulation.mat')
% load('e:\data\publications\2012\vikram_fact_des\code\codeOkPerm\ParamOptFactDes.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

visualisation = 1;
d = cell2mat(Dstat);
X = XX;
X2 = XX2;
Y = YY;

if visualisation,
    for i = 1:nFact
        figure
        y = d(i : nFact : end,:);
        surf(X,Y,y,'EdgeColor','none')
        hold on
        surf(X2,Y,y,'EdgeColor','none')
        xlabel('log_2(fold ratio)')
        ylabel('-log_1_0(p-value)')
        zlabel('SSQ')
        camlight left; lighting phong
        if simulation,
            title(['Factor ',num2str(i)])
        else
            title(['Factor ', factorLabel{i}])
        end
        axis tight
        colormap jet
        set(gca,'FontSize', FontSize)
        xl = get(gca,'title');
        set(xl,'FontWeight','normal');
        
        % log10(p-value) for each factor
        figure
        pValueMatrix = zeros(nRatio, nPValue);
        for j = 1:nRatio,
            for k = 1:nPValue,
                pValueMatrix(j,k)=PvalueTable{j,k}(i);
            end
        end
        pValueMatrix(find(pValueMatrix==0))=0.001;
        pValueMatrix = -1*log10(pValueMatrix);
        surf(X,Y,pValueMatrix,'EdgeColor','none')
        hold on
        surf(X2,Y,pValueMatrix,'EdgeColor','none')
        xlabel('log_2(fold ratio)')
        ylabel('-log_1_0(p-value)')
        zlabel('-log_1_0(p-value)')
        camlight left; lighting phong
        if simulation,
            title(['Factor ',num2str(i)])
        else
            title(['Factor ', factorLabel{i}])
        end
        axis tight
        colormap jet
        set(gca,'FontSize', FontSize)
        xl = get(gca,'title');
        set(xl,'FontWeight','normal');
    end

    if simulation,
        figure
        surf(X,Y,precision,'EdgeColor','none')
        hold on
        xlabel('log_2(fold ratio)')
        ylabel('-log_1_0(p-value)')
        zlabel('Precision')
        title('Precision')
        surf(X2,Y,precision,'EdgeColor','none')
        camlight left; lighting phong
        axis tight
        colormap jet
        set(gca,'FontSize', FontSize)
        xl = get(gca,'title');
        set(xl,'FontWeight','normal');
        
        figure
        surf(X,Y,recall,'EdgeColor','none')
        hold on
        xlabel('log_2(fold ratio)')
        ylabel('-log_1_0(p-value)')
        zlabel('Recall')
        title('Recall (sensitivity)')
        surf(X2,Y,recall,'EdgeColor','none')
        camlight left; lighting phong
        axis tight
        colormap jet
        set(gca,'FontSize', FontSize)
        xl = get(gca,'title');
        set(xl,'FontWeight','normal');
        
        figure
        surf(X,Y,gScore,'EdgeColor','none')
        hold on
        xlabel('log_2(fold ratio)')
        ylabel('-log_1_0(p-value)')
        zlabel('g-score')
        title('g-score')
        surf(X2,Y,gScore,'EdgeColor','none')
        camlight left; lighting phong
        axis tight
        colormap jet
        set(gca,'FontSize', FontSize)
        xl = get(gca,'title');
        set(xl,'FontWeight','normal');
        
        figure
        surf(X,Y,fScore,'EdgeColor','none')
        hold on
        xlabel('log_2(fold ratio)')
        ylabel('-log_1_0(p-value)')
        zlabel('f-score')
        title('f-score')
        surf(X2,Y,fScore,'EdgeColor','none')
        camlight left; lighting phong
        axis tight
        colormap jet
        set(gca,'FontSize', FontSize)
        xl = get(gca,'title');
        set(xl,'FontWeight','normal');
    end
    
    figure
    surf(X,Y,log10(sleectedPeaks),'EdgeColor','none')
    hold on
    xlabel('log_2(fold ratio)')
    ylabel('-log_1_0(p-value)')
    zlabel('log_1_0(number of selected peaks)')
    title('log_1_0 of number of selected variables')
    surf(X2,Y,log10(sleectedPeaks),'EdgeColor','none')
    camlight left; lighting phong
    axis tight
    colormap jet
    set(gca,'FontSize', FontSize)
    xl = get(gca,'title');
    set(xl,'FontWeight','normal');
end

%save all figures
figHandles = get(0,'Children');

mkdir('figureSaveDir');
for i = 1:length(figHandles),
    saveas(figHandles(i),['figureSaveDir\' num2str(i)]);
end