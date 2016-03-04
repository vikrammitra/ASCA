function [p_factor, p_interaction] = permutations(X, F, facs, interactions, PVth,FCth,center, n_perm, peaksAffected,  visualisation, sSize, factorLabel)

% Does permutation test on the factors of a crossed experimental design
%
% Input: 
% Data matrix X           rows: subject; columens: variables
% Design matrix F         rows: level per factor for each factor
% facs                    row vector with factors that must be included;
%                         the factors must correspond to the columns in F
% Interactions            array with per row interacting factors
% PVth                    therhold for t-test p-value for Volcano filtering
% FCth                    therhold for fold change for Volcano filtering
% center                  1: center 2: auto scale
% n_perm                  number of permutations
% peaksAffected           binary vetor to signal peaks affected by signiifcant factors
% visualisation           flag used for visualisation. 0 no plots, 1 visualise plots
% sSize                   number of factors
% factorLabel             cells contaning the labels of factors
% 
% Output:
% p_fcator                p-values vector for factors
% p_interaction           p-values vector for interactions
%
% GZ, 23 January 2012. Matlab 7.11.0
% Version 1.0
% Version 1.1  24 January 2012: added permutation test for interaction
% Version 1.2  1 Feb 2012: fixed estimate of interaction significance DJV
% Version 1.21 19 July 2012: fixed some comments
% GZ, 23 January, 2013
% Version 2.0
% - factors to evaluate can be chosen
% - improved output
%
% Copyright: Gooitzen Zwanenburg
% This code is available under APACHE Licence 2.0
% http://www.apache.org/licenses/LICENSE-2.0.html

tic
    [indx, XVolOrig] = selectionVolcano(X,F, PVth, FCth, peaksAffected, visualisation, factorLabel);
    size_data = size(XVolOrig);                              % size of the data matrix
    XVolOrig = dataCentralisation(XVolOrig, center);
    
    n_levels            = max(F);                               % number of levels / factor
    n_interactions      = size(interactions,1);                 % number of interactions
    factors             = facs;                                 % factors to be evaluated
    % determine if default factors are needed (default is all factors)
    if isempty(factors)
        n_factors = size(F,2);
        factors = 1 : n_factors;
        %disp('All factors are used');
    else
        n_factors       = size(factors,2);       % number of factors
        %disp('Factors used are:')
        %disp(factors)
    end
    if isempty(interactions)
        %disp('No interactions are used')
    else
        %disp('Interactions used:')
        %disp(interactions)
    end
    
    X_level_means       = cell(n_factors,1);                    % level_means per factor  
    X_permuted          = cell(n_factors,1);                    % permuted data per factor
    X_interaction_means = cell(n_interactions);                 % cell_means X_int - X_a - X_b
    ssq                 = zeros(n_perm + 1,n_factors);          % within sum of squares data matrix
    ssq_interaction     = zeros(n_perm + 1,n_interactions);     % within sum of squares data matrix
    p_factor            = zeros(1, n_factors);                  % p-values factors
    p_interaction       = zeros(1, n_interactions);             % p-values interactions
  
    % Make structure with unchanging 'variables'
    Fixed.F              = F;
    Fixed.interactions   = interactions; 
    Fixed.n_factors      = n_factors;
    Fixed.n_levels       = n_levels;
    Fixed.n_interactions = n_interactions;

% Collect level means for the factors indicated in the model 
    for factor = factors 
        X_level_means{factor} = level_means(XVolOrig, Fixed, factor);
    
        % sum of squares of level averages of measured data
        ssq(1, factor) = sum(sum( (X_level_means{factor}).^2));
    end
    
    for i = 1 : n_interactions                                            
        factor_1 = interactions(i,1);
        factor_2 = interactions(i,2);
        
        % Check that the factors of the interaction are calculated
        if (isempty(find(factors == factor_1, 1)) || isempty(find(factors == factor_2, 1)))
%             %disp('Warning, factor matrices for interaction not calculated');
%             %disp('Rerun with appropriate vector facs');
        end
        
        % X_residuals = X - Xa - Xb
        X_residuals = XVolOrig - X_level_means{factor_1} - X_level_means{factor_2};
        % Calculate cell means
       
        X_interaction_means{n_interactions} = zeros(size(XVolOrig));
        for level_factor_1 = 1 : n_levels(factor_1)          % levels for first factor             
            for level_factor_2 = 1 : n_levels(factor_2)      % levels for second factor  
                tmp = zeros(size_data(1),1);
                found = find((F(:,factor_2) == level_factor_2) & ...  % find rows 
                             (F(:,factor_1) == level_factor_1));
                m = mean(X_residuals(found,:));                         % average over cell   
                tmp(found) = 1;
                X_interaction_means{i} = X_interaction_means{i} + tmp*m;
            end
        end
        ssq_interaction(1, i) = sum(sum( (X_interaction_means{i}).^2));
    end   
    
    
    % Interactions p-values are calculated through permutation of X - Xa -
    % Xb
    % Do the permutations (do this, for now, for two factors)
    for j = 1 : n_perm

        % Permutations of experimental factors
        % Restrict permutations to the factor considered
        for factor = factors
%             permuted = zeros(size_data(1),1);
        
            % only rows for which the levels of "the other factors" do not
            % change may be permuted, i.e. the rows for which the elements for
            % each column in the F-matrix, except for column 'factor', are the
            % same can be permuted.
        
            F_reduced = Fixed.F;
            F_reduced(:,factor) = F_reduced(randperm(length(F_reduced(:,factor))),factor);
            FixedPerm = Fixed;
            FixedPerm.F = F_reduced;
            [indx, XVolPerm] = selectionVolcano(X,F_reduced, PVth, FCth, peaksAffected, visualisation, factorLabel);
            size_data = size(XVolPerm);
            XVolPerm = dataCentralisation(XVolPerm, center);
            X_permuted{factor} = XVolPerm;

            % Level means for each factor 
            X_level_means{factor} = level_means(X_permuted{factor}, FixedPerm, factor);
   
            ssq(1 + j, factor) = sum(sum( (X_level_means{factor}).^2));
        end
        
        % Permutations for interactions. 
        for i = 1 : n_interactions               
            factor_1 = interactions(i,1);
            factor_2 = interactions(i,2);
            F_reduced = Fixed.F;
            F_reduced(:,interactions(i,1)) = F_reduced(randperm(length(F_reduced(:,interactions(i,1)))),interactions(i,1));
            F_reduced(:,interactions(i,2)) = F_reduced(randperm(length(F_reduced(:,interactions(i,2)))),interactions(i,2));
            FixedPerm = Fixed;
            FixedPerm.F = F_reduced;
            [indx, XVolPerm] = selectionVolcano(X,F_reduced, PVth, FCth, peaksAffected, visualisation, factorLabel);
            size_data = size(X_r);
            X_r = dataCentralisation(X_r, center);

            % Calculate level means experimental factors for permuted data 
            for factor = [factor_1 factor_2]
                X_level_means{factor} = level_means(X_r, FixedPerm, factor);
            end
            
            X_residuals = X_r - X_level_means{factor_1} - X_level_means{factor_2};
            X_interaction_means{n_interactions} = zeros(size(X));
            for level_factor_1 = 1 : n_levels(factor_1)                   % levels for first factor             
                for level_factor_2 = 1 : n_levels(factor_2)               % levels for second factor  
                    tmp = zeros(size_data(1),1);
                    found = find((F(:,factor_2) == level_factor_2) & ...  % find rows 
                                 (F(:,factor_1) == level_factor_1));
                    m = mean(X_residuals(found,:));                       % average over cell   
                    tmp(found) = 1;
                    X_interaction_means{i} = X_interaction_means{i} + tmp*m;
                end
            end
            ssq_interaction(1 + j, i) = sum(sum( (X_interaction_means{i}).^2));
        end           
%         disp(j)
    end
    
    % Calculate p-values
    % how many ssq's are larger than data ssq?
    for factor = factors
        p_factor(factor) = size(find(ssq(2:n_perm + 1, factor) ...
                        >= ssq(1, factor)),1)/(n_perm);
    end
    for interaction = 1 : n_interactions
        p_interaction(interaction) = size(find(ssq_interaction(2:n_perm + 1, interaction) ...
                                  >= ssq_interaction(1, interaction)),1)/(n_perm);
    end
    
%      %disp('=====================================================');
%      %disp('=                     Results                       =')
%      %disp('=====================================================');
%      
%      fprintf('\n\n')
%      fprintf('Factor: ');
%      fprintf('%8d', factors);
%      fprintf('\nP-value:');
%      fprintf('%8.3f',p_factor(factors));
%      fprintf('\nSSQ%:', ssq);
%      fprintf('\n\n')
 toc   
end            % function

function M = level_means(Y, Input, factor)
    % Calculate level means for factor in matrix Y
    size_Y = size(Y);
    M = zeros(size_Y);
    for level = 1 : Input.n_levels(factor)
        tmp = zeros(size_Y(1),1);
        found = find(Input.F(:,factor) == level);  % find rows that belong to level
        m = mean(Y(found,:));                      % calculate level mean 
        tmp(found) = 1;                            % flag the rows found
        M = M + tmp*m;
    end
end

function XMean = dataCentralisation(X,center)
    size_data = size(X);                              % size of the data matrix
    if center == 1
        XMean = (X - ones(size_data(1),1)*mean(X));      % Centre
    elseif center == 2
        XMean = (X - ones(size_data(1),1)*mean(X))./...  % Standardize
        (ones(size_data(1),1)*std(X));
    end
end

function [indx, XSelected] = selectionVolcano(X,F, PVth, FCth, peaksAffected, visualisation, factorLabel)
    % center/standardize the data
    [indx] = volcano_plot(X, F, PVth, FCth, peaksAffected, 0, 0, visualisation, factorLabel);
    if length(indx)>0,
        XSelected = X(:,indx);
    else
        XSelected = X;
    end
    XSelected = log10(XSelected);
    % chek if data contains 0. If yes warning and relplace -Inf by 0.
    if ~(isempty(find(XSelected == -Inf))),
        disp('-Inf value! for log)')
        XSelected(find(XSelected == -Inf))=0;
    end
end