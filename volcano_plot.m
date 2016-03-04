function [imp_dat] = volcano_plot(X,F,alpha,th_t,peaksAffected,logR,maxPRatio,visualisation, legendString)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Volcano plot and peak filtering
%
% Authors: Peter Horvatovich, Vikram Mitra and Gooitzen Zwanenburg
% Date: July 17, 2015
%
% Input variables:
% X:                quantitative data matrix
% F:                low/high level distribution of factors
% alpha:            -log10 p-value of 2 sample independent t-test threshold
% th_t:             log2 fold ratio threshold
% peaksAffected:    boolean vector where 1 mark affected peaks, if it is empty all peaks visualised with dot
% logR:             1 calcualte only ratio, other valie log2 ratio
% maxPRatio:        0 make selection, 1 does not make selection
% visualisation:    boolean 0 non-visualisation, 1 Volcano plot
% legendString:     cell with the name of factors used in Volcano plot
%
% Output variable:
% imp_dat:          index of selected peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nFact = size(F,2);
sSize = size(F,1);
FontSize = 14;

size_data = size(X);
X_level_means  = cell(nFact,1);
n_levels = size(unique(F),1); 
m = cell(n_levels,nFact);
fold_rat = zeros(size_data(2),nFact);
pvalues = zeros(size_data(2),nFact);
colorFact=[0,0,1;0,1,0;1,0,1;...
    1,0,0;0.26,0.85,1;1,0.85,0.11;...
    0.36,0.75,0.26];
if visualisation,
    figure
    hold on
    set(gca,'FontSize', FontSize)
    title(['Peak selection with -log_1_0(p-value): ', num2str(alpha), ' and log_2(fold ratio): ', num2str(th_t(1)), ' and ', num2str(th_t(2))])
    xlabel('log_2(fold change)')
    ylabel('-log_1_0(p-value)')
    legendHandles = zeros(nFact,1);
end
for factor = 1:nFact
    X_level_means{factor} = zeros(size_data);
    
    for level = 1 : n_levels
        Dummy = find(F(:,factor) == level);      % find rows that belong to level
        m{level,factor} = X(Dummy,:);            % level data
    end
    
    high = m{1,factor}';
    low = m{2,factor}';
    %For each factor get the level differences
    d = high - low;
    if logR == 1;
        ratio = (mean(high,2)./mean(low,2));
    else
        ratio = log2(mean(high,2)./mean(low,2));
    end
    if size(find(ratio ==-Inf|ratio ==Inf| isnan(ratio)),1)>0,
        ratio(find(ratio ==-Inf|ratio ==Inf| isnan(ratio)))=0;
        disp('Ratio problem');
    end
    fold_rat(:,factor) = ratio;
    [h,dummyP] = ttest2(high', low');
    pvalues(:,factor)=-1*log10(dummyP);
    
    %select the points that are within the thershold
    if ~maxPRatio,
        [imp_indx_t] = find(ratio < th_t(1) | ratio > th_t(2));
        [imp_indx_p] = find(pvalues(:,factor)> alpha);
        imp_indx{factor} = intersect(imp_indx_t,imp_indx_p);
        
        if visualisation,
            if length(peaksAffected)>0,
                legendHandles(factor) = scatter(fold_rat(find(peaksAffected(factor,:)==0),factor),pvalues(find(peaksAffected(factor,:)==0),factor),'.'...
                    ,'MarkerFaceColor',colorFact(factor,:),'MarkerEdgeColor',colorFact(factor,:),'SizeData',20);
                scatter(fold_rat(find(peaksAffected(factor,:)),factor),pvalues(find(peaksAffected(factor,:)),factor),'+'...
                    ,'MarkerFaceColor',colorFact(factor,:),'MarkerEdgeColor',colorFact(factor,:),'SizeData',46);
                scatter(fold_rat(imp_indx{factor},factor),pvalues(imp_indx{factor},factor),'o'...
                    ,'MarkerFaceColor',colorFact(factor,:),'MarkerEdgeColor',colorFact(factor,:),'MarkerFaceColor','none','SizeData',46);
            else
                legendHandles(factor) = scatter(fold_rat(:,factor),pvalues(:,factor),'b.','MarkerEdgeColor',colorFact(factor,:),'MarkerFaceColor','none','SizeData',20);
                scatter(fold_rat(imp_indx{factor},factor),pvalues(imp_indx{factor},factor),'o'...
                    ,'MarkerFaceColor',colorFact(factor,:),'MarkerEdgeColor',colorFact(factor,:),'MarkerFaceColor','none','SizeData',46);
            end
        end
    end
end

if visualisation,
    if length(legendString)>0,
        legend(legendHandles, legendString);
    end
end
if maxPRatio,
    imp_dat = [max(max(abs(pvalues))), max(max(abs(fold_rat)))];
else
    imp_dat=[];
    for i=1:nFact,
        imp_dat = [imp_dat; imp_indx{i}];
    end
    imp_dat=unique(sort(imp_dat));
end