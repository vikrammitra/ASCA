% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ASCA analysis of simulated data
% %
% % Authors: Peter Horvatovich, Vikram Mitra and Gooitzen Zwanenburg
% % Date: July 17, 2015
% %
% % set mainDir variable to exact location of the ASCA main directory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% clc
% clear
% warning off
% % change the directory according to the exact coce location
% mainDir = 'e:\data\publications\2012\vikram_fact_des\code\codeOkPerm\';
% % mainDir = 'd:\users\peter\munka\publications\2012\vikram_fact_des\code\codeOk\';
% cd(mainDir)
% nRun = 100; % number of repetation
% factorLabel = {'\bfFactor 1 ','Factor 2 ','\bfFactor 3 ','Factor 4 ','\bfFactor 5 ','Factor 6 ','Factor 7 '};
% 
% f = importdata('FMatrix.txt'); % experiment design file (the same as the original fractional factorial design)
% f([4 11 14],:) = [];
% f = f([2,7,5,14,1,15,3,6,12,10,11,8,13,16,4,9],:);
% nParam = 2559; % number of variables
% sSize = 16; % number of samples
% nFact = 7; % number of factors
% const = 6.043; % constant for simulation and for replace 0 in the real data
% sd = 0.5391; % standard deviation for simulation
% beta=1;
% lowM = [-1.5,0,2,0,-3,0,0]; % mean for lower level
% highM = [1.5,0,-2,0,3,0,0]; % mean for higher level
% percentagePeaks = [0.05, 0, 0.035, 0, 0.05, 0, 0]; % % of factor affected peaks
% peaksNAffected=round(repmat(nParam, nFact,1).*percentagePeaks')'; % vector contaning the number of factor influenced peaks
% visualisation = 0; % use 1 to visualize Volcano plot and 0 to not.
% FontSize = 14; % Fontsize in the plots
% 
% % test ASCA with Fixed pvalue and Fold change thresholds
% th1 = [0; -log2(2); -2;-2];
% th2 = -th1;
% th = [th1,th2];
% 
% pth = [0; -log10(0.05); -log10(0.05); 2];
% Dstat = cell(size(pth,1),nRun);
% PvalueTable = cell(size(pth,1),nRun);
% peaksAffectedNos = zeros(size(pth,1),nRun);
% 
% peaksAffected=cell(size(pth,1),nRun); % matrix containing the index of peaks, which are influenced by the factor (with 1)
% Recall = zeros(size(pth,1),nRun);
% precision = zeros(size(pth,1),nRun);
% selectedPeaks = zeros(size(pth,1),nRun);
% contingencyTable = cell(size(pth,1),nRun);
% specificity = zeros(size(pth,1),nRun);
% gScore = zeros(size(pth,1),nRun);
% fScore = zeros(size(pth,1),nRun);
% IDX = cell(size(pth,1),nRun);
% 
% % setting parallel execution
% nPool = 7; % number of threads
% poolobj = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(poolobj)
%     parpool('local',nPool);
% end
% 
% for k = 1:size(pth,1),
%     parfor i = 1:nRun,
%         disp(['Iteration: ' num2str(k) ' and ' num2str(i)]);
%         R_main = zeros(sSize ,nParam);
%         peaksAffected{k,i} = zeros(nFact,nParam);
%         
%         for nn=1:nFact,
%             Dummy=randperm(nParam);
%             peaksAffected{k,i}(nn, Dummy(1,1:peaksNAffected(nn)))=ones(1,peaksNAffected(nn));
%         end
%         peaksUniqueAffected=max(peaksAffected{k,i});
%         tindx = find(peaksUniqueAffected);
%         nindx = find(peaksUniqueAffected==0);
%         
%         for mn = 1:nFact,
%             % peaks affceted
%             high = find(f(:,mn) == 2);
%             low = find(f(:,mn) == 1);
%             if peaksNAffected(mn)>0,
%                 R_main(low,find(peaksAffected{k,i}(mn,:)))=R_main(low,find(peaksAffected{k,i}(mn,:)))+normrnd(lowM(mn),sd,[size(low,1),peaksNAffected(mn)]);
%                 R_main(high,find(peaksAffected{k,i}(mn,:)))= R_main(high,find(peaksAffected{k,i}(mn,:)))+normrnd(highM(mn),sd,[size(high,1),peaksNAffected(mn)]);
%             end
%             R_main(:,find(peaksAffected{k,i}(mn,:)==0))= R_main(:,find(peaksAffected{k,i}(mn,:)==0))+normrnd(0,sd,sSize,[size(find(peaksAffected{k,i}(mn,:)==0),2)]);
%         end
%         
%         peaksAffectedNos(k,i) = length(tindx);
%         R_main = R_main/sqrt(nFact);
%         R_main = R_main + const;
%         R_main_10p = exp(R_main);
%         
%         if find(R_main == 0),
%             disp('0 in the matrix');
%         end;
%         
%         % Test ASCA with variying fold change
%         [XASCA,indx,p_dat] = ASCAwithPerm(R_main_10p,f,sSize, [1:nFact], pth(k,:),th(k,:),peaksAffected{k,i},visualisation, factorLabel);
% %         [indx] = volcano_plot(R_main_10p,f,pth(k,:),th(k,:),peaksAffected{k,i},0,0,visualisation, factorLabel);
%         selectedPeaks(k,i)=size(indx,1);
%         indxNotSel = 1:nParam;
%         indxNotSel(indx) = [];
%         ASCA{k,i} = XASCA;
%         D = XASCA.effects.ssq_factors';
% %         if length(indx)>0,
% %             XASCA=R_main_10p(:,indx);
% %         else
% %             XASCA=R_main_10p;
% %         end
%         if ~(isempty(XASCA)),
% %             % ASCA on the entire data
% %             log_n = log10(XASCA);
% %             if ~(isempty(find(log_n == -Inf))),
% %                 disp('-Inf value! for log)')
% %                 log_n(find(log_n == -Inf))=0;
% %             end
% %             if size(log_n,2)>=sSize,
% %                 ASCA{k,i} = asca(log_n,f,[],[],1,visualisation);
% %                 D = ASCA{k,i}.effects.ssq_factors';
% %                 for Fpval = 1:size(f,2)
% %                     PvalueFactors = permutations(log_n,f,[Fpval],[],1,10000);
% %                     p_dat{Fpval,1} = PvalueFactors(Fpval);
% %                 end
% %             else
% %                 D=zeros(nFact,1);
% %                 p_dat = mat2cell(ones(nFact,1))
% %             end
%             Dstat{k,i} = D(:,1);
%             PvalueTable{k,i} = p_dat;
%             IDX{k,i} = indx;
%             % tindx true indeces and nindx are false indeces
%             tp = length(intersect(tindx,indx));
%             fp = length(intersect(nindx,indx));
%             tn = length(intersect(nindx,indxNotSel));
%             fn = length(intersect(tindx,indxNotSel));
%             contingencyTable{k,i}=[tp,fp;fn,tn];
%             if (tp+fp+tn+fn)~=nParam,
%                 disp('Problem with contigency table!')
%             end
% 
%             % Calculate specificity and sensitivity
%             recall(k,i) = tp/(tp+fn);
%             precision(k,i) = tp/(tp+fp);
%             specificity(k,i) = tn/(tn+fp);
%             gScore(k,i)= sqrt(specificity(k,i)*recall(k,i));
%             fScore(k,i)=((1+beta^2)*(precision(k,i)*recall(k,i)))/((beta^2)*precision(k,i)+recall(k,i));
%         else
%             Dstat{k,i} = zeros(nFact,1);
%             recall(k,i) = 0;
%             precision(k,i) = 0;
%             IDX{k,i} = [];
%             PvalueTable{k,i} = ones(nFact,1);
%         end
%     end
% end
% delete(poolobj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualisation
% 
% This code part is used to load saved simulated data and make visualisation
clc
clear
%sload('e:\data\publications\2012\vikram_fact_des\code\codeOkPerm\Simulation100.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


visualisation = 1;
PvalueTableMatrix=cell(4,1);
if visualisation,
    for i = 1:size(pth,1),
%         PvalueTableMatrix=cell(4,1);
        SSQ = [];
        for j = 1:nRun,
            PvalueTableMatrix{i}=[PvalueTableMatrix{i};PvalueTable{i,j}];
            SSQ = [SSQ; Dstat{i,j}'];
        end
        PvalueTableAverage(i,:) = -log10(mean(PvalueTableMatrix{i}));
        SSQAverage(i,:) = mean(SSQ);
        PvalueTableStdev(i,:) = std(PvalueTableMatrix{i});
        SSQStdev(i,:) = std(SSQ);
        
        %SSQ bar plot
        figure
        hold on
        set(gca,'FontSize', FontSize)
        barwitherr(SSQStdev(i,:), SSQAverage(i,:), 0.5, 'FaceColor', [0 0.553 0.914], 'LineStyle', 'none');
        ylabel('SSQ%')
        set(gca, 'XTickLabel',factorLabel, 'XTick',1:numel(factorLabel), 'LineWidth', 1.0)
        rotateXLabels(gca, 90)
        dummyAxis = axis;
        dummyAxis = [dummyAxis(1)+0.5, dummyAxis(2)-0.5, dummyAxis(3) 1.1*max(SSQAverage(i,:))];
        axis(dummyAxis)
        if pth(i)>0,
            title(['SSQ%: peaks selected with ' num2str(pth(i)) ' log_1_0(p-value) and ±' num2str(th(i,2)) ' log_2(fold ratio)'], 'FontSize', FontSize-3)
        else
            title('SSQ%: using all peaks', 'FontSize', FontSize-2)
        end
        xl = get(gca,'title');
        set(xl,'FontWeight','normal');

        %p-value bar plot
        figure
        hold on
        set(gca,'FontSize', FontSize)
        barwitherr(PvalueTableStdev(i,:), PvalueTableAverage(i,:), 0.5, 'FaceColor', [1 0.25 0.25], 'LineStyle', 'none');
        ylabel('-log_1_0(p-value)')
        set(gca, 'XTickLabel',factorLabel, 'XTick',1:numel(factorLabel))
        rotateXLabels(gca, 90)
        dummyAxis = axis;
        dummyAxis = [dummyAxis(1)+0.5, dummyAxis(2)-0.5, 0 1.1*max(PvalueTableAverage(i,:))];
        axis(dummyAxis)
        line([0, dummyAxis(2)],[-log10(0.05), -log10(0.05)], 'LineStyle', '--', 'LineWidth', 0.75)
        if pth(i)>0,
            title(['SSQ -log_1_0(p-value): peaks selected with ' num2str(pth(i)) ' log_1_0(p-value) and ±' num2str(th(i,2)) ' log_2(fold ratio)'], 'FontSize', FontSize-3)
        else
            title('SSQ -log_1_0(p-value): using all peaks', 'FontSize', FontSize-2)
        end
        xl = get(gca,'title');
        set(xl,'FontWeight','normal');
    end
end

% peak selected in all repetition
loadingsHigh = zeros(nFact,nRun);
loadingsLow = zeros(nFact,nRun);
NaffcetedPeaks = zeros(nFact,nRun);
NNonAffcetedPeaks = zeros(nFact,nRun);
expN = 2;
for i = 1:nRun
    for k = 1:nFact
        [dummyBool,dummyIdx] = ismember(IDX{expN,i},find(peaksAffected{expN,i}(k,:)),'R2012a');
        NaffcetedPeaks(k,i) = length(find(dummyIdx));
        dummyLoadings = -1*abs(ASCA{expN,i}.factors.loadings{k}(:,1));
        [rankASCA,TIEADJASCA] = tiedrank(dummyLoadings);
        loadingsHigh(k,i) = mean(rankASCA(dummyBool));
        dummyIdxLoadings = 1:length(dummyBool);
        dummyIdxLoadings(dummyBool) = [];
        NNonAffcetedPeaks(k,i) = length(dummyIdxLoadings);
        loadingsLow(k,i) = mean(rankASCA(dummyIdxLoadings));
    end
end

%save all figures
figHandles = get(0,'Children');
figureSaveDir = [mainDir 'Simulation100\'];
for i = 1:length(figHandles),
    saveas(figHandles(i),[figureSaveDir 'Sim100Figure_' num2str(i)]);
end