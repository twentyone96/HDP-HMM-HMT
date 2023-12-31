% GenHypo lists M most possible hypotheses, calculates their 
%   probabilities and form the cell array cellHypoNew.
%   
% Syntax: 
%   cellHypoNew = GenHypo(cellHypoSeed, cellTargSeed, M,...
%       Z, H, R, pd, densNew, densClt)
%   
% In: 
%   cellHypoSeed - a cell array containing father (seed) hypotheses and 
%       their probabilities. Each cell contains both a hypo and its 
%       log-prob. A hypo is a set of (meas, targ) assignment.  
%   cellTargSeed - a cell array containing description of targets. 
%       Each cell contains information of a full collection of targets
%       derived according to the corresponding hypo. One target was
%       described with {idx, lifePoint, X, P}. 
%   M - the number of hypoes needed to generate using Murty's algorithm
%   Z - a dMeas*numMeas matrix containing meas on this step
%   H - measurement matrix
%   R - measurement covariance matrix
%   pd - target detection probability
%   densNew - density of new target in the surveillance space
%   densClt - density of clutter in the surveillance space  
%       
% Out: 
%   cellHypoNew - a cell array containing newly generated hypoes and 
%       their probabilities

function cellHypoNew = GenHypo(cellHypoSeed, cellTargSeed, M,...
    Z, H, R, pd, densNew, densClt)

global maxMes

len = size(cellHypoSeed, 2);
%cellHypoNew = [];
asso = size(maxMes,1);
prob = zeros(1,1);
hypo = {asso, prob};
cellHypoNew = cell(1,M*len);
% for i = 1 : M*len
%     cellHypoNew{i} = {hypo};
% end

for i = 1 : len
    %% 航迹初始化阶段
    if isempty(cellTargSeed{i})
        cellTargSeed{i} = {[]};
    end
    targ = cellTargSeed{i}; % {nTarg * {idx, lifePoint, X, P}}
                            % in accord with cellHypoSeed{i}
    if ~isempty(targ{1})
%		if(targ{3}{1}{1} ~=  [0,0,0,0,0,0]')
        maxTargIdx = max(cellfun(@(v) v{1}, targ));
    else
        maxTargIdx = 0;
    end
    
    nMeas = size(Z, 2);
    
    % form the log-prob matrix 
    ProbMat = GenProbMat(maxTargIdx, targ, ...
        Z, H, R, pd, densNew, densClt);
    
    % use Murty's algorithm to find out M minimum assignments
    cellHypoTmp = Murty(ProbMat, M); %列是所有观测，行是目标/新目标/杂波，每一个观测只能属于一列（一种情况）
%     cellHypoTmp = Murty(ProbMat, 1);
%     cellHypoTmp = {cellHypoTmp{1},cellHypoTmp{1},cellHypoTmp{1},cellHypoTmp{1}};
    hypoSeed = cellHypoSeed{i};
    if ~isempty(hypoSeed{1})
        p1 = hypoSeed{2}; % log-prob of this hypoSeed
        nT = sum(hypoSeed{1} > 0); % number of targets deteced in hypoSeed
    else
        p1 = 0;
        nT = 0;
    end
    p2 = -nT * log(1-pd); 
    
    for j = 1 : M
        hypo = cellHypoTmp{j};
        asso = hypo{1};
        
        % set indicators of clutter to 0 in assignment matrix
        % see mat1, mat2 and mat3 in GenProbMat
        asso(asso > maxTargIdx + nMeas) = 0;
        
        % make indicators to new targets more compact
        idxNew = find(asso > maxTargIdx);
        asso(idxNew) = maxTargIdx+1:maxTargIdx+length(idxNew); 

        % go on with probability calculation
        prob = hypo{2};
        prob = prob + p1 + p2; 
        prob = round(prob,3);
        
        hypo{1} = asso;
        hypo{2} = prob;
        cellHypoTmp{j} = hypo;
        cellHypoNew{(i-1)*M + j} = cellHypoTmp{j};
    end
%    cellHypoNew = [cellHypoNew, cellHypoTmp];
end

