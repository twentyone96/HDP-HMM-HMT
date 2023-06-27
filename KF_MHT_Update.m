%% 这个版本在预测步和更新步都用了复合估计
% KF_MHT_Update generates target states according to multiple hypotheses
% 
% Syntax:
%   cellTargNew = KF_MHT_Update(cellTargSeed, cellHypoNew,...
%       M, Z, H, R, maxLifePoint)
%   
% In: 
%   cellTargSeed - a cell array containing description of targets. 
%       Each cell contains information of a full collection of targets
%       derived according to the corresponding hypo. One target was
%       described with {idx, lifePoint, X, P}. 
%   cellHypoNew - a cell array containing newly generated hypoes and 
%       their probabilities. 
%       Note: size(cellHypoNew, 2) == M * size(cellTargSeed, 2)
%   M - number of hypoes generated using Murty's algorithm
%   Z - measurements at this time step
%   H - measurement matrix
%   R - measurement covariance
%   maxLifePoint - maximum life point: each target maintains a "life 
%       point", at each time step, if this target is not detected, then
%       subtract 1 from its life point; else if its life point is less
%       than maxLifePoint, then add 1 to its life point.
%   
% Out: 
%   cellTargNew - a cell array containing newly generated targets info
%       that is one-one correspondent to cellHypoNew. Each cell in 
%       cellHypoNew leads to one cell in cellTargNew which gives a full
%       set of state estimations for all the targets.  
%
% Description: 
%  cellHypo{i} is {asso_i, prob_i}
%  cellTarg{i} is {nTarg_i*{targInfo}}, 
%      and targInfo{i} is {idx_i, lifePoint_i, X_i, P_i}

function cellTargNew = KF_MHT_Update(cellTargSeed, cellHypoNew,...
    M, Z, H, R, maxLifePoint, chooseIdx, t)
global MaxStateNum
global dim

rear = 0;
cellTargNew = cell(1, size(cellHypoNew, 2));
idx = zeros(1,1);
liftPoint = zeros(1,1);
targ = cell(1,19);
targ{1} = zeros(1,6);
targ{2} = eye(6,6);
targ{3} = zeros(1,1,'int8');
targ{4} = zeros(1,1,'int8');
targ{5} = zeros(2,1);
targ{6} = zeros(MaxStateNum,MaxStateNum);
targ{7} = zeros(1,1,'int8');
targ{8} = zeros(1,1+MaxStateNum,'int8');
targ{9} = zeros(1,1);
targ{10} = zeros(MaxStateNum,MaxStateNum,'int8'); 
targ{11} = zeros(1,1);
targ{12} = zeros(1,1);
targ{13} = zeros(MaxStateNum,1);
targ{14} = zeros(1,1);
targ{15} = zeros(dim,1);
targ{16} = zeros(1,1);
targ{17} = zeros(dim,dim);
targ{18} = zeros(MaxStateNum,1);
targ{19} = zeros(MaxStateNum,4);
% targs = cell(1,100);
% targs{1:100} = {targ};

%aTarg = {idx, liftPoint, targs};
%cellTargNew{1:size(cellHypoNew, 2)} = {aTarg};

updateIdx = 0;
%% X(1) P(2) currentState(3) L(4) cellCtr(5)...
%% cellTrans(6) alpha(7) beta(8) gamma(9) M(10)...
%% auxPhi(11) auxTau(12) auxZeta(13)
%% kappa(14) vartheta(15)(dim=u的维度，这里选2位） nu(16) Delta(17)（2*2矩阵）
%% --这一行都是NIW的参数
%% NIW(0.001, 0, 50, I_u(2*2)) 
%% cellStates(18) 每个状态保持的数量

for i = 1 : size(cellTargSeed, 2)
    head = rear + 1;
    rear = rear + M;
    for j = head : rear
        updateIdx = updateIdx + 1;
        %% 只有在范围内的才计算，计算量可是100个粒子
         aCase = {[]};
         if ismember(updateIdx, chooseIdx) == 1
%             aCase = {[]};
%         else
            % use cellHypoNew{j} to update cellTargSeed{i}, thus form
            % cellTargNew{j}
            asso = cellHypoNew{j}{1};   % 取关联数据
            if isempty(cellTargSeed{i})
                cellTargSeed{i} = {[]};
%								cellTargSeed{i} = aTarg;
            end
            
            if ~isempty(cellTargSeed{i}{1})
%						if(cellTargSeed{i}{3}{1}{1} ~=  [0,0,0,0,0,0]')
                nExistedTarg = size(cellTargSeed{i}, 2); 
                maxTargIdx = max(cellfun(@(v) v{1}, cellTargSeed{i}));
            else
                nExistedTarg = 0; 
                maxTargIdx = 0;
            end
            nNewTarg = sum(asso > maxTargIdx);

            aCase = cell(1, nExistedTarg+nNewTarg);
            % 1. deal with existed (in last step) targets 
            for k = 1 : nExistedTarg % for each target
                aTarg = cellTargSeed{i}{k}; % one target
                idx = aTarg{1}; % the index of aTarg
                lifePoint = aTarg{2}; 
 %               X = aTarg{3};
 %               P = aTarg{4};

                if lifePoint == 0 % a disappeared target
%                    aTarg{3} = {aTarg{3}{1}};
                    aCase{k} = aTarg;
                    continue; % just pass it
                end
                flg = find(asso == idx);
                if isempty(flg) % there is no meas asso with aTarg
                    lifePoint = lifePoint - 1;
                else
                    aMeas = Z(:, flg); % the meas asso with aTarg
                    aTarg = UpdateExitTarg(aTarg, aMeas, t);
                    % Kalman filter update stage
%                     innov = aMeas - H*X;
%                     S = R + H*P*H';
%                     G = P*H'/S;
%                     X = X + G*innov;
%                     P = P - G*S*G';

                    if lifePoint < maxLifePoint
                        lifePoint = lifePoint + 1;
                    end
                end

                aTarg{2} = lifePoint;
%                 aTarg{3} = X;
%                 aTarg{4} = P;

                aCase{k} = aTarg;
            end
            % 2. deal with newly observed tentative targets
            for k = 1 : nNewTarg
                idx = maxTargIdx + k;
                flg = find(asso == idx);
                aMeas = Z(:, flg); 

                % initialize a new target
                aTarg = cell(1, 3);
                aTarg{1} = idx;
                aTarg{2} = 3;%1;% 3; % lifePoint  或者设置为1
                newTarg = NewTarg(aMeas);
                aTarg{3} = newTarg;
%                 aTarg{3} = [aMeas(1), 0, aMeas(2), 0]'; % X
%                 aTarg{4} = diag([1 1 1 1]); % P

                aCase{idx} = aTarg;
            end
        end
        cellTargNew{j} = aCase;
    end
end





