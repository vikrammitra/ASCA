function loopParamOpt(th, pth, R_main_10p, f, nFact, visualisation, factorLabel, nParam, kkk, sSize,simulation, FontSize, nRatio, nPValue, mainDir)
% function to perform parameter optimisation
% main reason to place this part of code in function is to allow forlop
% parallelisation
% 
% Parameters:
% th            fold change threshold
% pth           p-value threshold
% R_main_10p    Original data matrix
% f             levels distribution according to factors
% nFact         number of factors
% visualisation flag for visualisation (0 no visualisation (off); 1 visuliastion on)
% factorLabel   name of factors
% nParam        number of variables
% kkk           interation number
% sSize         number of samples
% simulation    flag to define is simulated data or real data is used. 1 means simulation, 0 real data
% FontSize      font size of text that are 
% nRatio        number of ratio threshold to be tested
% nPValue       number of p-value threshold to be tested
% mainDir       main working directory

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
            disp(['Iteration: ' num2str(i) ' and ' num2str(j) ' and ' num2str(kkk)]);
            if ~(isempty(XASCA)),
                Dstat{i,j} = D(:,1);
                PvalueTable{i,j} = p_dat;
            else
                Dstat{i,j} = zeros(nFact,1);
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
    figureSaveDir = [mainDir 'ParameterOptFactdesEmpty\'];
    for i = 1:length(figHandles),
        saveas(figHandles(i),[figureSaveDir 'ParamOptFactDesFigure_' num2str(i) '_' num2str(kkk)],'meta');
    end
    save([mainDir 'ParameterOptFactdesEmpty\emptySetSimulation_' num2str(kkk) '.mat'])
    close all
end