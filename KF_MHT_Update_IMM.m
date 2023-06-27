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

function cellTargNew = KF_MHT_Update_IMM(cellTargSeed, cellHypoNew,...
    M, Z, H, R, maxLifePoint, T)
global dim
global IMMNUM

rear = 0;
cellTargNew = cell(1, size(cellHypoNew, 2));

%aTarg = {idx, liftPoint, targs};
%cellTargNew{1:size(cellHypoNew, 2)} = {aTarg};

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
%        updateIdx = updateIdx + 1;
        %% 只有在范围内的才计算，计算量可是100个粒子
%         if ismember(updateIdx, chooseIdx) == 1
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
				        cellIMM = aTarg{3};
				%       X = zeros(3*dim*IMMNUM, 1);
				%       P = zeros(dim*3*IMMNUM,dim*3);
				%       modeP = zeros(1,IMMNUM);
				%       coeffi = zeros(1,IMMNUM);
				%       zinv = zeros(dim,IMMNUM);
				%       zcov = zeros(dim,IMMNUM);
				        
					      X = cellIMM{1};
					      P = cellIMM{2};
					      coeffi = cellIMM{3};
					      modeP = cellIMM{4};

 %              X = aTarg{3};
 %              P = aTarg{4};

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
%                    aTarg = UpdateExitTarg(aTarg, aMeas, t);
                    for ii=1:IMMNUM
                        if X(3*dim*(ii-1)+2) == 0
                            X(3*dim*(ii-1)+2) = (aMeas(1) - X(3*dim*(ii-1)+1))/T;
                            X(3*dim*(ii-1)+5) = (aMeas(2) - X(3*dim*(ii-1)+4))/T;
                            X(3*dim*(ii-1)+1) = aMeas(1);
                            X(3*dim*(ii-1)+4) = aMeas(2);
                            X(3*dim*(ii-1)+3) = aMeas(3);
                            X(3*dim*(ii-1)+6) = aMeas(3);
                        else
                            % Kalman filter update stage
                              innov = aMeas(1:2) - H*X(3*dim*(ii-1)+1:3*dim*ii);
                              S = R + H*P(dim*3*(ii-1)+1:dim*3*ii,1:3*dim)*H.';
                              G = P(dim*3*(ii-1)+1:dim*3*ii,1:3*dim)*H.'/S;
                              X(3*dim*(ii-1)+1:3*dim*ii) = X(3*dim*(ii-1)+1:3*dim*ii)...
                                + G*innov;
                              P(dim*3*(ii-1)+1:dim*3*ii,1:3*dim) = ...
                                P(dim*3*(ii-1)+1:dim*3*ii,1:3*dim) - G*S*G.';
                        end
                    end
                    if lifePoint < maxLifePoint
                        lifePoint = lifePoint + 1;
                    end
                end
                aTarg{2} = lifePoint;
%                 aTarg{3} = X;
%                 aTarg{4} = P;
                cellIMM{1} = X;
                cellIMM{2} = P;
                cellIMM{3} = coeffi;
                cellIMM{4} = modeP;

                aTarg{3} = cellIMM;

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
                aTarg{2} = 3; % lifePoint  或者设置为1
                X = [aMeas(1), 0, 0,aMeas(2), 0,0].';
                X = [X;X;X;X];
                P = diag([1 1 1 1 1 1]);
                P = [P;P;P;P];
                modeP = [1;0.00;0.00;0.00];
                coeffi = [1;0.00;0.00;0.00];
                zinv = zeros(dim,IMMNUM);
                zcov = zeros(dim,IMMNUM);
                cellIMM{1} = X;
                cellIMM{2} = P;
                cellIMM{3} = coeffi;
                cellIMM{4} = modeP;

                aTarg{3} = cellIMM;
%                newTarg = NewTarg(aMeas);
%                aTarg{3} = newTarg;
%                 aTarg{3} = [aMeas(1), 0, aMeas(2), 0]'; % X
%                 aTarg{4} = diag([1 1 1 1]); % P

                aCase{idx} = aTarg;
            end
%        end
        cellTargNew{j} = aCase;
    end
end






