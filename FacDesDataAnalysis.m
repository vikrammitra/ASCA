% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ASCA analysis of experimental design data
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
% load('factorLabel.mat')
% factorLabel{2} = ['\bf' factorLabel{2}];
% factorLabel{5} = ['\it' factorLabel{5}];
% factorLabel{6} = ['\bf' factorLabel{6}];
% for i = 1:length(factorLabel),
%     factorLabel{i} = [factorLabel{i} ' '];
% end
% nRun = 100; % number of repetation
% 
% f = importdata('FMatrix.txt'); % experiment design file (the same as the original fractional factorial design)
% f([4 11 14],:) = [];
% f = f([2,7,5,14,1,15,3,6,12,10,11,8,13,16,4,9],:);
% sSize = 16; % number of factors
% nFact = 7;
% const = 6.043; % constant for simulation and for replace 0 in the real data
% sd = 0.5391; % standard deviation for simulation
% FontSize = 14;% Fontsize in the plots
% visualisation = 0; % use 1 to visualize Volcano plot and 0 to not.
% FontSize = 12; % Fontsize in the plots
% mzFilter = 0; %was 400
% rtFilterLow = 45; %was 50
% rtFilterHigh = 135;%was 13
% 
% % test ASCA with Fixed pvalue and Fold change thresholds
% th1 = [0; -log2(2); -2;-2];
% th2 = -th1;
% th = [th1,th2];
% 
% pth = [0; -log10(0.05); -log10(0.05); 2];
% Dstat = cell(size(pth,1),nRun);
% PvalueTable = cell(size(pth,1),nRun);
% 
% data = importdata('Outstem_mzRadius=0.3_TRadius=1_Fraction=0.5_mzStart=100_rtStart=65_mzEnd=1500_rtEnd=135.mpks');
% % removing data points < 400da Mz and < 50Mins and >130Mins
% data_FS = data.data;
% [filtdata_mz indices] = find(data_FS(:,2) > mzFilter);
% filtdata_mz = data_FS(filtdata_mz,:);
% [filtdata_rt indices] = find(filtdata_mz(:,3) > rtFilterLow & filtdata_mz(:,3) <rtFilterHigh);
% filtdata_rt = filtdata_mz(filtdata_rt,:);
% dat_FS = filtdata_rt(:,14:end)';
% nParam = size(dat_FS,2); % number of variables
% selectedNPeaks = zeros(size(pth,1),nRun);
% selectedPeaks = cell(size(pth,1),nRun);
% 
% for k =1:size(pth,1),
%     for i = 1:nRun,
%         disp(['Iteration: ' num2str(k) ' and ' num2str(i)]);
%         % resampling and replacing missing intensity values
%         NonZeroMat = dat_FS;
%         for ii = 1:size(NonZeroMat,2)
%             zeroIn = find(NonZeroMat(:,ii) == 0);
%             nZeros = length(zeroIn);
%             R = lognrnd(const,sd,[size(NonZeroMat,1),1]);
%             NonZeroMat(zeroIn,ii) = R(zeroIn);
%         end
%         R_main_10p = NonZeroMat;
%         
%         %Test ASCA with variying fold change
%         [XASCA,indx,p_dat] = ASCAwithPerm(R_main_10p,f,sSize, [1:nFact], pth(k,:),th(k,:),[],visualisation, factorLabel);
%         selectedNPeaks(k,i) = size(indx,1);
%         selectedPeaks{k,i} = indx;
%         indxNotSel = 1:nParam;
%         indxNotSel(indx) = [];
%         ASCA{k,i} = XASCA;
%         D = XASCA.effects.ssq_factors';
%         ASCALoadings{k,i} = XASCA.factors.loadings;
% %         [indx] = volcano_plot(R_main_10p,f,pth(k,:),th(k,:),[],0,0,visualisation, factorLabel);
% %         selectedNPeaks(k,i) = size(indx,1);
% %         selectedPeaks{k,i} = indx;
% %         indxNotSel=1:nParam;
% %         indxNotSel(indx)=[];
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
% %                 ASCA = asca(log_n,f,[],[],1,visualisation);
% %                 ASCALoadings{k,i} = ASCA.factors.loadings;
% %                 D = ASCA.effects.ssq_factors';
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
%         else
%             Dstat{k,i} = zeros(nFact,1);
%             PvalueTable{k,i} = ones(nFact,1);
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualisation
% 
% This code part is used to load saved simulated data and make visualisation
% clc
% clear
% load('e:\data\publications\2012\vikram_fact_des\code\codeOkPerm\FactDes100.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

visualisation = 0;
if visualisation,
    for i = 1:size(pth,1),
        PvalueTableMatrix=cell(size(pth));
        PvalueTableAverage=[];
        PvalueTableStdev=[];
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
        ylabel('SSQ')
        set(gca, 'XTickLabel',factorLabel, 'XTick',1:numel(factorLabel), 'LineWidth', 1.0)
        rotateXLabels(gca, 90)
        dummyAxis = axis;
        dummyAxis = [dummyAxis(1)+0.5, dummyAxis(2)-0.5, dummyAxis(3) 1.1*max(SSQAverage(i,:))];
        axis(dummyAxis)
        if pth(i)>0,
            title(['SSQ: peaks selected with ' num2str(pth(i)) ' log_1_0(p-value) and ±' num2str(th(i,2)) ' log_2(fold ratio)'],'FontSize',FontSize-2)
        else
            title('SSQ: using all peaks','FontSize',FontSize-2)
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
            title(['-log_1_0(p-value): peaks selected with ' num2str(pth(i)) ' log_1_0(p-value) and ±' num2str(th(i,2)) ' log_2(fold ratio)'],'FontSize',FontSize-2)
        else
            title('-log_1_0(p-value): using all peaks','FontSize',FontSize-2)
        end
        xl = get(gca,'title');
        set(xl,'FontWeight','normal');
    end
end

% peak selected in all repetition
nExp = 1;
intersectSelectedPeaks = selectedPeaks{nExp,1};
rankFact = cell(nFact,1);
TIEADJ = cell(nFact,1);
FacDesDataIdx = [];
for i = 1:nRun,
    intersectSelectedPeaks = intersect(intersectSelectedPeaks,selectedPeaks{nExp,i});
end
for i = 1:nRun
    [dummyC, dummyia, dummyib] = intersect(intersectSelectedPeaks,selectedPeaks{nExp,i});
    for k = 1:nFact
        loadingsASCA{k,i} = abs(ASCALoadings{nExp,i}{k}(dummyib,1));
        [rankFact{k}(i,:),TIEADJ{k}(i,:)] = tiedrank(-1*loadingsASCA{k,i});
    end
end
for k = 1:nFact
    rankProductColumns = prod(rankFact{k},1);
    [sortedRanksProduct{k},sortedRankIdx{k}] = sort(rankProductColumns, 'ascend');
    FacDesDataIdx(k,:) = intersectSelectedPeaks(sortedRankIdx{k});
end

%Annotation selected peaks
factSel = [2,5,6];
nPeakSelected = 10;
load([mainDir 'FinalAnnotationVectorOk.mat'])
proteinAnnotation  = readtable([mainDir 'upload\annotatedProteinOkOk.xlsx'],'FileType','spreadsheet');
FacDesDataAnnotation = zeros(size(FacDesDataIdx));
FacDesDataPeptideAnnotation = cell(nFact,nPeakSelected);
FacDesDataAnnotationLabel = cell(nFact,nPeakSelected);
FacDesDataPeptideProteinAnnotationLabel = cell(nFact,nPeakSelected);
mzTable = zeros(nFact,nPeakSelected);
rtTable = zeros(nFact,nPeakSelected);
for k = 1:nFact
    for kk = 1:length(FacDesDataIdx),
        dummyAnnotation = matchedPeaksIdx(find(FacDesDataIdx(k,kk)==matchedPeaksIdx(:,1)),2);
        if length(dummyAnnotation)>0,
            FacDesDataAnnotation(k,kk)=dummyAnnotation;
        end
    end
    for kk = 1:nPeakSelected,
        if FacDesDataAnnotation(k,kk),
            dummyString = annotationTable.proteinID{FacDesDataAnnotation(k,kk)};
            dummyStringPeptide = annotationTable.peptideSeq{FacDesDataAnnotation(k,kk)};
%             dummyString = strsplit(dummyString,'|');
            FacDesDataAnnotationLabel{k,kk} = dummyString;
            FacDesDataPeptideAnnotation{k,kk} = dummyStringPeptide;
            mzTable(k,kk) = annotationTable.mz(FacDesDataAnnotation(k,kk));
            rtTable(k,kk) = annotationTable.rt(FacDesDataAnnotation(k,kk));
            dummyString3 = strcat(table2array(proteinAnnotation(k,kk)), ['\n' dummyStringPeptide]);
            FacDesDataPeptideProteinAnnotationLabel{k,kk} = dummyString3{1};
        else
            FacDesDataAnnotationLabel{k,kk} = 'not annotated';
            FacDesDataPeptideAnnotation{k,kk} = 'not annotated';
            FacDesDataPeptideProteinAnnotationLabel{k,kk} = 'not annotated';
        end
    end
end

%get mz and rt and index from the original raw data
discPeaks = cell(length(factSel),1);zeros(length(factSel),nPeakSelected);
loadingsSelected = cell(length(factSel),1);
for i = 1:length(factSel),
    discPeaks{i} = filtdata_rt(FacDesDataIdx(factSel(i),1:nPeakSelected),[1:3 9]);
    loadingsSelected{i} = zeros(nRun,nPeakSelected);
    for k = 1:nRun,
        [dummyC, dummyia, dummyib] =  intersect(FacDesDataIdx(factSel(i),1:nPeakSelected),selectedPeaks{nExp,k},'stable');
        loadingsSelected{i}(k,:) = ASCALoadings{nExp,k}{factSel(i)}(dummyib,1);
    end
end

visualisation = 0;
%create bar plots for loadings
for i = 1:length(factSel),
    figure
    dummyMean = mean(abs(loadingsSelected{i}));
    dummySd = std(abs(loadingsSelected{i}));
    hold on
    set(gca,'FontSize', FontSize)
    barwitherr(dummySd, dummyMean, 0.5, 'FaceColor', [0 0.553 0.914], 'LineStyle', 'none');
    axisValue = axis;
    axisValue([1 2]) = [0.5 10.5];
    axis(axisValue);
%     set(gca,'Xtick',1:nPeakSelected,'XTickLabel',FacDesDataPeptideProteinAnnotationLabel(factSel(i),:),'FontSize', FontSize-2)
    set(gca,'Xtick',1:nPeakSelected,'XTickLabel',table2array(proteinAnnotation(factSel(i),:)),'FontSize', FontSize-2)
    rotateXLabels(gca, 90);
%     format_ticks(gca,FacDesDataPeptideProteinAnnotationLabel(factSel(i),:),[],1:nPeakSelected,[],90,90,0.1)
    ylabel('Average loadings')
    title(factorLabel{factSel(i)})
end

%create box plots for differential analysis
% subplotV=[2,5];
% subplotFil=[1,2,3,4,5,6,7,8,9,10];
BoxPlotlocation = 'e:\data\publications\2012\vikram_fact_des\code\codeOkPerm\Boxplots\';
for i = 1:length(factSel),
    for k = 1:nPeakSelected,
        h=figure
        set(h, 'DefaultTextFontSize', FontSize+8);
        dummyIdx = FacDesDataIdx(factSel(i),k);
        dataHigh = dat_FS(find(f(:,factSel(i))==1),dummyIdx);
        dataLow = dat_FS(find(f(:,factSel(i))==2),dummyIdx);
        boxplot([dataHigh,dataLow],'labels',{'Low','High'});
        hold on
        scatter(repmat(1,1,length(dataHigh)),dataHigh,175,'.','MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.8 0])
        scatter(repmat(2,1,length(dataLow)),dataLow,175,'.','MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.8 0])
        set(gca,'FontSize', FontSize+6)
%         title([factorLabel{factSel(i)},' ',num2str(filtdata_rt(dummyIdx,1)),' index, ',num2str(filtdata_rt(dummyIdx,2),'%8.3f'),' m/z, ',...
%             num2str(filtdata_rt(dummyIdx,3),'%6.3f'),' rt ', num2str(filtdata_rt(dummyIdx,4),'%2i'),' peaks'])
        saveas(h, [BoxPlotlocation 'Boxplot_fact_' factorLabel{factSel(i)}(4:end-1) '_' num2str(k) '_peak.fig'],'fig');
        saveas(h, [BoxPlotlocation 'Boxplot_fact_' factorLabel{factSel(i)}(4:end-1) '_' num2str(k) '_peak.emf'],'meta');
        close(h)
    end
end

% % This code should be used with TAPP to create input file for export EICs
% % use TIC_visualisation.m script to visualise EICs obtained with TAPP
% 
% % create EIC file for TAPP
% TAPPDir = [mainDir 'TicVisualisation\'];
% inputFileName = [TAPPDir 'alingedSpectraList.txt'];
% outputFileNameMesh = [TAPPDir 'inputTIC_16FactDes.txt'];
% outputFileNameWarp = [TAPPDir 'inputTIC_16FactDes_warp.txt'];
% fileID = fopen(inputFileName, 'r');
% fileList = textscan(fileID,'%s', 'Delimiter','\n');
% fclose(fileID);
% cleanOrigList = cell(size(fileList{1}));
% tmapList = cell(size(fileList{1}));
% for i = 1:size(fileList{1},1),
%    dummy = strsplit(fileList{1}{i},' ');
%    dummy2 = strsplit(dummy{1},'.');
%    tmapList{i} = dummy2{1};
%    dummy3 = strsplit(dummy2{1},'_');
%    cleanOrigList{i} = dummy3{2};
% end
% clear dummy
% 
% tickness = 0.25;
% startTime = 45;
% endTime = 135;
% 
% fileID = fopen(outputFileNameMesh,'w');
% for i=1:size(fileList{1},1),
%     fprintf(fileID,'%s\n',['File ' TAPPDir cleanOrigList{i} '.mesh']);
%     for k = 1:length(factSel),
%         for kk = 1:nPeakSelected,
%             dummyIdx = FacDesDataIdx(factSel(k),kk);
%             fprintf(fileID,'%s\n',[num2str(filtdata_rt(dummyIdx,2)),' ', num2str(tickness), ' ', num2str(startTime), ' ', num2str(endTime)]);
%         end
%     end
% end
% fclose(fileID);
% 
% fileID = fopen(outputFileNameWarp,'w');
% for i=1:size(fileList{1},1),
%     fprintf(fileID,'%s\n',[TAPPDir tmapList{i} '.tmap ']);
% end
% fclose(fileID);
% 
% %save all figures
% figHandles = get(0,'Children');
% figureSaveDir = [mainDir 'FactDes100\'];
% for i = 1:length(figHandles),
%     saveas(figHandles(i),[figureSaveDir 'Sim100Figure_' num2str(i)]);
% end