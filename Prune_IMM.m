% Prune cuts less possible hypoes and corresponding state estimations 
% from cellHypo and cellTarg as well as take record of good estimations. 
% 
% Syntax: 
%   [cellEstm, cellHypo, cellTarg] = ...
%       Prune(cellHypo, cellTarg, M, N, t, LP)
% 
% In:
%   cellEstm - a cell array containing the state estimation for each 
%       target. The structure of each cell is {idx, startTime, matX}.
%       The startTime is the time step on which the target is confirmed.
%   cellHypo - a cell array of size 1*(1+M+M^2+...+M^N). The structure 
%       of each cell in the cell array is {asso, prob}, where asso is a
%       column vector denoting the assignment of a set of measurements; 
%       prob is a scalar denoting narual log of likelihood of this hypo.
%   cellTarg - a cell array of the same size of cellHypo. Each cell in 
%       the cell array is called oneCase which gives a set of estimations 
%       for all targets according to the correspoding hypo in cellHypo. 
%       The structure of each cell in oneCase is {idx, lifePoint, X, P} 
%       describing one target. 
%   M - number of hypoes generated for one seed hypo
%   N - scan depth
%   t - current time step. Note MHT introduces time delay, so at time 
%       step t we can only confirm things happened at step t-N+1. 
%   LP - the confirmed hypo may contains some tentative targets that are
%       not included in previous estimations (i.e. input cellEstm). If
%       the lifePoints of such targets are not less than LP, then they
%       are confirmed to be new targets. 
%
% Out: 
%   cellEstm - updated cellEstm
%   cellHypo - cellHypo after pruning
%   cellTarg - cellTarg after pruning  
  
function [cellEstm, cellHypo, cellTarg] = Prune_IMM...
    (cellEstm, cellHypo, cellTarg, M, N, t, LP)

% use the latest probabilities as the standard to judge
head = 2;
for i = 1 : N-1
    head = head + M^i;
end
rear = head + M^N - 1;      % 第一次运行，元素22~85
arrayProb = cellfun(@(v) v{2}, cellHypo(head : rear));  % 提取最后一层的prob
x = find(arrayProb == min(arrayProb)); % -log(prob), so find minimum 概率小，-log(prob)越大，所以找最小的概率最大
chooseBranch = ceil(x/M^(N-1)); % 每一个root有4个分枝，每一个分枝末端是M^(N-1)即16个节点
idx = 1;
chooseIdx = [];
for i = 1 : N
    chooseIdx = [chooseIdx, ...
        idx+(chooseBranch-1)*M^(i-1)+1 : idx+chooseBranch*M^(i-1)]; % 2~5代表最基本的4个假设分枝
    idx = idx + M^i;
end

% pick up the right branch
cellHypo = cellHypo(chooseIdx);
cellTarg = cellTarg(chooseIdx);

% update cellEstm
confirmedCase = cellTarg{1}; % the root of the hypo tree
for i = 1 : size(confirmedCase, 2)
    aTarg = confirmedCase{i};
    idx = aTarg{1};
    lifePoint = aTarg{2};
    %% 所有的目标都存储在
    if lifePoint == 0
       continue;
    else
        cellIMM = aTarg{3};
        flg = [];
        if ~isempty(cellEstm)
            flg = find(cellfun(@(v) v{1}, cellEstm) == idx);
        end
        if isempty(flg) % maybe a new target appears    目标起始为0，这个目标就是新目标
            if lifePoint >= LP
                newTarg = {idx, t-N+1, cellIMM};
%                if t - N + 1 > 4
%                    X = {X{1}};
%                    newTarg = {idx, t-N+1, X};
%                    cellEstm = [cellEstm, {newTarg}];   
%                else
%                    newTarg = {idx, t-N+1, X};
                    cellEstm = [cellEstm, {newTarg}];
%                end
            end
        else % an existed target
            thisTarg = cellEstm{flg};
            thisTarg{3} = cellIMM;
            cellEstm{flg} = thisTarg;
        end
    end
end





