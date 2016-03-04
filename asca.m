function ASCA = asca(X, F, facs, interactions, center, visualisation, display)

% ASCA does a principal component analysis on the level averages of each
% experimental factor in a designed experiment with balanced data.
% Interactions between two factors can also be calculated
% Nomenclature: 
% Factor : experimental factor such as time or dose
% Level  : treatment belonging to a factor, for example three time points
%          correspond to three levels of the factor time, or two doses
%          correspond to two levels of the factor dose.
% Cell   : combination of two factors/levels         

% The program takes a data matrix X and a matrix F that describes the
% experimental setup as input.

% Input
% data matrix X  :   each row is a measurement, each column a variable
% design matrix F:   columns correspond to factors rows to levels. This
%                    matrix must be constructed by the researcher.
% facs           :   row vector with factors that must be included
%                    the factors must correspond to the columns in F
% interactions   :   vector with per row the factors for which interactions
%                    are to be calculated
% center         :   center == 0 use raw data; center == 1: center data; 
%                    center == 2 standardize data
% visualisation  :   visualisation == 0 not plots are provided;
%                    visualisation == 1 figures plot are provided
% display        :   display == 0 no display (faster execution)
%                :   display == 1 display (slower execution)
%
% Output
% structure ASCA that contains scores, loadings singular values and
% projections of the
% factors and interactions
%
% GZ, 17 January 2012. Matlab 7.11.0
% Version 1.0
% GZ, 22 May 2012. Matlab 7.11.0
% Version 1.1
% - Added calculation of percentage explained
% - Added plotting of the interactions and their residuals
% GZ, 18 June, 2012
%Version 1.2
% - Allowed for cells with single subject. This caused a problem when 
%   calculating the interaction means. Put in an 'if' statement to catch
%   single subject cells.
% GZ, 2 October, 2012
%Version 1.3
% - Added percentage 'explained' for overall mean, factors, interactions
% and residuals.
% GZ, 23 January, 2013
% Version 2.0
% - factors to evaluate can be chosen
% - improved output
% - added some input checking
%
% Copyright: Gooitzen Zwanenburg
% This code is available under APACHE Licence 2.0
% http://www.apache.org/licenses/LICENSE-2.0.html

    clear ASCA;
    size_data           = size(X);                   % size of the data matrix 
    n_levels            = max(F);                    % number of levels / factor
    n_interactions      = size(interactions,1);      % number of interactions
    factors             = facs;                      % factors to be evaluated
    
    % determine if default factors are needed (default is all factors)
    if isempty(factors)
        n_factors = size(F,2);
        factors = 1 : n_factors;
        if display,
            disp('All factors are used');
        end
    else
        n_factors       = size(factors,2);       % number of factors
        if display
            disp('Factors used are:')
            disp(factors)
        end
    end
    if isempty(interactions)
        if display
            disp('No interactions are used')
        end
    else
        if display
            disp('Interactions used:')
            disp(interactions)
        end
    end
  
    X_raw               = X;                         % Input data matrix
    X_level_means       = cell(n_factors,1);         % level_means per factor  
    SSQ_factors         = zeros(1,max(factors),1);      % sum of squares for factors
    X_interaction_means = cell(n_interactions);      % cell_means for interactions 
    SSQ_interactions    = zeros(1,n_interactions);   % sum of squares for interactions
    
    % Structure with results
    ASCA.factors.scores               = cell(n_factors,1);
    ASCA.factors.loadings             = cell(n_factors,1);
    ASCA.factors.projected            = cell(n_factors,1);
    ASCA.factors.singular             = cell(n_factors,1);
    ASCA.factors.explained            = cell(n_factors,1);
    ASCA.interactions.scores          = cell(n_interactions);
    ASCA.interactions.loadings        = cell(n_interactions); 
    ASCA.interactions.singular        = cell(n_interactions); 
    ASCA.interactions.explained       = cell(n_interactions); 
    

    % center/standardize the data
    if center == 1
        Mean = ones(size_data(1),1)*mean(X_raw);        % Overall mean
        X = (X_raw - Mean);                             % Center
        SSQ_mean = sum(sum(Mean.^2));                   % SSQ overall means
        SSQ_X = sum(sum(X_raw.^2));                     % Sum of squares data matrix
    elseif center == 2
        Mean_std = ones(size_data(1),1)*mean(X)./...  
                  (ones(size_data(1),1)*std(X));
        X_std = X_raw./(ones(size_data(1),1)*std(X_raw));
        X = (X_std - Mean_std);                         % Standardize
        SSQ_mean = sum(sum(Mean_std.^2));               % SSQ overall means
        SSQ_X = sum(sum(X_std.^2)); 
    end
    X_residuals         = X;                            % initial matrix with residuals    
    ASCA.data           = X;
    ASCA.design.F       = F;
    ASCA.design.facs    = factors;

    % Collect level means for the factors indicated in the model 
    for factor = factors 
        X_level_means{factor} = zeros(size_data);
        for level = 1 : n_levels(factor)
            tmp = zeros(size_data(1),1);
            found = find(F(:,factor) == level);      % find rows that belong to level       
            m = mean(X(found,:));                    % calculate level mean 
            tmp(found) = 1;                          % flag the rows found
            X_level_means{factor} = X_level_means{factor} + tmp*m;
        end
        SSQ_factors(factor) = sum(sum(X_level_means{factor}.^2));
        X_residuals = X_residuals - X_level_means{factor};
    end
    
    % Interactions
    for i = 1 : n_interactions                                            
        factor_1 = interactions(i,1);
        factor_2 = interactions(i,2);
        
        % Check that the factors of the interaction are calculated
        if (isempty(find(factors == factor_1, 1)) || isempty(find(factors == factor_2, 1)))
            if display
                disp('Warning, factor matrices for interaction not calculated');
                disp('Rerun with appropriate vector facs');
            end
            return
        end
        
        X_interaction_means{i} = zeros(size_data);
        for level_factor_1 = 1 : n_levels(factor_1)          % levels for first factor             
            for level_factor_2 = 1 : n_levels(factor_2)      % levels for second factor 
                tmp = zeros(size_data(1),1);
                found = find((F(:,factor_2) == level_factor_2) & ...  % find rows 
                             (F(:,factor_1) == level_factor_1));
                if size(found,1) == 1                        % only one subject/cell
                    m = X(found,:);                          % average is row
                else
                    m = mean(X(found,:));                    % average over cell
                end
                tmp(found) = 1;
                X_interaction_means{i} = X_interaction_means{i} + tmp*m;
            end
        end
        X_interaction_means{i} = X_interaction_means{i} - ...
                                 X_level_means{factor_1} - X_level_means{factor_2};
        SSQ_interactions(i) = sum(sum(X_interaction_means{i}.^2));                     
        X_residuals = X_residuals - X_interaction_means{i};                     
    end   
     
    SSQ_residuals = sum(sum(X_residuals.^2));
    perc_effects = effect_explains(SSQ_X, SSQ_mean, SSQ_factors, ...
                                   SSQ_interactions, SSQ_residuals);
    ASCA.effects = perc_effects;
    ASCA.residuals = X_residuals;
    %Do PCA on level averages for each factor
    for factor = factors
        [t,s,p] = do_pca(X_level_means{factor});
        projected_data = (X_residuals*p);       % project residuals on loadings
        ASCA.factors.scores{factor}    = t;
        ASCA.factors.loadings{factor}  = p;
        ASCA.factors.singular{factor}  = s;
        ASCA.factors.projected{factor} = projected_data;
        ASCA.factors.explained{factor} = pc_explains(s);
    end

    %Do PCA on interactions
    for interaction = 1 : n_interactions
        [t,s,p] = do_pca(X_interaction_means{interaction});
        projected_data = (X_residuals*p);       % project residuals on loadings
        ASCA.interactions.scores{interaction}    = t;
        ASCA.interactions.loadings{interaction}  = p;
        ASCA.interactions.singular{interaction}  = s;
        ASCA.interactions.projected{interaction} = projected_data;
        ASCA.interactions.explained{interaction} = pc_explains(s);
    end

    if visualisation,
        plot_asca(ASCA);
        plot_interactions(ASCA, n_interactions);
    end
    
    if display
        disp('Percentage each effect contributes to the total sum of squares')
        disp('Overall means')
        disp(perc_effects.ssq_mean);
        disp('Factors')
        disp(factors);
        disp(perc_effects.ssq_factors(perc_effects.ssq_factors >0))
        disp('Interactions')
        disp(perc_effects.ssq_interactions)
        disp('Residuals')
        disp(perc_effects.ssq_residuals)
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [scores, singular_values, loadings] = do_pca(D)
        % This function does the work: do PCA through singular value
        % decomposition on input matrix. Returns scores, loadings and
        % singular values.
        [u, sv, v] = svd(D);
        scores = u*sv;
        singular_values = diag(sv);
        loadings = v;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    function plot_asca(ASCA)

    % Function to plot the scores, loadings and projections for each 
    % factor.
    
    % Input
    % strucuture ASCA

    % Output
    % none
    
        f_F              = ASCA.design.F; 
        f_n_levels       = max(f_F);                    % number of levels / factor
        f_n_factors      = ASCA.design.facs;            % factors
        legend_labels    = cell(max(f_n_levels),1);
        
    % Plot results for experimental factors
        for f_factor = f_n_factors
            loadings       = ASCA.factors.loadings{f_factor};
            scores         = ASCA.factors.scores{f_factor};
            projections    = ASCA.factors.projected{f_factor};
            plot_residuals = scores + projections;
            figure;
            bar(loadings(:,1:2));
            legend('First loading', 'Second loading');
            sload = ['Loadings factor: ' int2str(f_factor)];
            title(sload);
            figure;
            hold;
            j = 1;                         % Count legends
            clear legend_labels;
            for ii = 1 : f_n_levels(f_factor)
                f_found = find(f_F(:,f_factor) == ii);
                if ii == 1
                    sym = '*';
                    co = 'b';
                    clr = [sym co];
                elseif ii == 2
                    sym = '*';
                    co = 'r';
                    clr = [sym co];
                elseif ii == 3
                    sym = '*';
                    co = 'g';
                    clr = [sym co];
                elseif ii == 4
                    sym = '*';
                    co = 'k';
                    clr = [sym co];
                elseif ii == 5
                    sym = '*';
                    co = 'y';
                    clr = [sym co];
                end
                legend_labels{j} = ['Score level ' int2str(ii)];
                legend_labels{j + 1} = ['Level ' int2str(ii)];
                j = j + 2;
                plot(scores(f_found,1), scores(f_found,2), clr, 'LineWidth', 3 );
                plot(plot_residuals(f_found,1), plot_residuals(f_found,2), clr);
                stitle = [ 'Score plot factor: ' int2str(f_factor)];
                title(stitle);
                legend(legend_labels);
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_interactions(ASCA, i_interactions)
        
    % Function to plot the scores, loadings and projections for each 
    % interaction.
    
    % Input
    % strucuture ASCA
    % number of interactions

    % Output
    % none
    
    % Plot results for experimental factors
        for i_inter = 1 : i_interactions
            loadings       = ASCA.interactions.loadings{i_inter};
            scores         = ASCA.interactions.scores{i_inter};
            projections    = ASCA.interactions.projected{i_inter};
            figure;
            bar(loadings(:,1:2));
            legend('First loading', 'Second loading');
            sload = ['Loadings interaction: ' int2str(i_inter)];
            title(sload);
            figure;        
            hold
            plot(scores(:,1), scores(:,2), '*r', 'LineWidth', 3 );
            plot(projections(:,1), projections(:,2), 'ob');
            % We want the text shifted with respect to the points
            x_axis_size = xlim;
            y_axis_size = ylim;
            x_dist = abs(x_axis_size(1) - x_axis_size(2));
            y_dist = abs(y_axis_size(1) - y_axis_size(2));
            % Shift text 2% up and to the right
            text(scores(:,1) + 0.02*x_dist, scores(:,2) + 0.02*y_dist, ...
                int2str(ASCA.design.F), 'FontSize', 8 );
         
            
            stitle = [ 'Score plot interaction: ' int2str(i_inter)];
            title(stitle);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function perc_explained = pc_explains(sv)
 
    % Function to calculate the percentage of variance explained by each PC
    
    % Input
    % vector with singular values 
    
    % Output
    % vector with percentages explained variance
    
        sv_squared     = sv.^2;     
        total_variance = sum(sv_squared);
        perc_explained = (sv_squared/total_variance)*100;
    end
    
    function perc_explained_effect = effect_explains(ssq_X, ssq_mean, ssq_factors, ...
                                     ssq_interactions, ssq_residuals)
        
        % Function to calculate the percentage of variance explained by
        % each effect.
        
        % Input
        % sums of squares of the data matrix, the factors and the
        % interactions
        
        % Output
        % vector with percentages explained variance for each effect.
        
        perc_explained_effect.ssq_mean = 100*(ssq_mean/ssq_X);
        perc_explained_effect.ssq_factors = 100*(ssq_factors./ssq_X);
        perc_explained_effect.ssq_interactions = 100*(ssq_interactions./ssq_X);
        perc_explained_effect.ssq_residuals = 100*(ssq_residuals/ssq_X);
    end
         
end    