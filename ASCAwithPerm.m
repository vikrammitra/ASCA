function [ASCARes,indx,p_dat, SSQ, SSQInt] = ASCAwithPerm(X, f, sSize, factors, PVth, FCth, peaksAffected, visualisation, factorLabel)
% ASCA calculation with outer loop permutation for claculation of p-value
%
% Parameters:
% input:
% X                 full data Matrix.
% f                 factors levels distribution matrix.
% sSize             number of samples
% factors           vector contaning the factors to be studied.
% PVth              p-value threshold for Volcano selection.
% FCth              Fold Change for Volcano selection.
% peaksAffected     matrix with 1 to indicated peak affected by factors. Use only in simulation.
% visualisation     boolean flag for visualisation of Volcano plot. 1 visualisation on, 0 visualisation off.
% factorLabel       facors' label. Used in Volcano plot visualisation.
% 
% output:
% XASCA             ASCA results structure
% indx              index of peaks selected after Volcano filtering
% p_dat             factors p-value
% SSQ               relative SSQ (normalised to the total SSQ) of factors
% SSQInt            relative SSQ (normalised to the total SSQ) of
% interactions
% 
% by Peter Horvatovich, 29 January 2016

% ASCA using original fcator levels
[indx] = volcano_plot(X, f, PVth, FCth, peaksAffected, 0, 0, visualisation, factorLabel);
nPerm = 1000;
disp(['P-value threshold:' num2str(PVth)])
disp(['Fold change threshold:' num2str(FCth)])

if length(indx)>0,
    XASCA = X(:,indx);
else
    XASCA = X;
end
if ~(isempty(XASCA)),
    % ASCA on the entire log10 data
    logXASCA = log10(XASCA);
    % chek if data contains 0. If yes warning and relplace -Inf by 0.
    if ~(isempty(find(logXASCA == -Inf))),
        disp('-Inf value! for log)')
        logXASCA(find(logXASCA == -Inf))=0;
    end
    if size(logXASCA,2)>=sSize,
        ASCARes = asca(logXASCA,f,[],[],1,visualisation, 1);
        % This function is using ASCA function to calcualte permutation
        % based p-value. It is a generic and can be used for any design,
        % but the function is slow.
        % p_dat = permutationPValue(X,f,factors,nPerm,PVth,FCth, peaksAffected, visualisation, sSize, factorLabel, ASCARes.effects.ssq_factors);

        % this function is used to calculate permutation based p-value
        % using 2 level designs. It is faster than the generic solution but
        % can be only used when number of levels is limited to 2 for all
        % factors.
        [p_dat, p_interaction] = permutations(X,f,factors,[],PVth,FCth,1,nPerm,peaksAffected, visualisation, sSize, factorLabel);
    else
        ASCARes.effects.ssq_factors = zeros(length(factors),1);
        p_dat = ones(length(factors),1);
    end
end

function [p_factor] = permutationPValue(X, f, factors, nPerm, PVth, FCth, peaksAffected, visualisation, sSize, factorLabel, SSQOrig)



SSQ = -ones(nPerm,length(factors));
for i = 1:nPerm,
    parfor ii = factors,
        fPerm = f;
        fPerm(:,ii) = fPerm(randperm(length(fPerm(:,ii))),ii);
        [indx] = volcano_plot(X,fPerm,PVth,FCth,peaksAffected,0,0,visualisation, factorLabel);
        if length(indx)>0,
            XASCA = X(:,indx);
        else
            XASCA = X;
        end
        % ASCA on the entire log10 data
        logXASCA = log10(XASCA);
        % chek if data contains 0. If yes warning and relplace -Inf by 0.
        if ~(isempty(find(logXASCA == -Inf))),
            disp('-Inf value! for log)')
            logXASCA(find(logXASCA == -Inf))=0;
        end
        if size(logXASCA,2)>=sSize,
            XASCA = asca(logXASCA,fPerm,[],[],1,visualisation, 0);
            SSQ(i,ii) = XASCA.effects.ssq_factors(ii);
        else
            SSQ(i,ii) = 0;
        end
    end
    disp(i)
end
for i = factors,
    p_factor(i) = length(find(SSQ(i,1:nPerm)>= SSQOrig(i)))/(nPerm);
end