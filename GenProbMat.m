% GenProbMat forms the cost matrix, or log-prob matrix here
% 
% Syntax:
%     ProbMat = GenProbMat(maxTargIdx, targ, ...
%         Z, H, R, pd, densNew, densClt)
%     
% In: 
%   maxTargIdx - the maximum index of target that had been associtated 
%       with a measurement
%   targ - a cell array {nTarg * {idx, lifePoint, X, P}}. Each cell in 
%       targ is a cell describing one specific target.
%   Z - a dMeas*nMeas matrix containing meas on this step
%   H - measurement matrix
%   R - measurement covariance matrix
%   pd - probability of detection
%   densNew - density of new target in the surveillance space
%   densClt - density of clutter in the surveillance space  
%   
% Out: 
%   ProbMat - a  matrix containing the log probabilities of each 
%       assocition. 
      
function ProbMat = GenProbMat(maxTargIdx, targ, ...
    Z, H, R, pd, densNew, densClt)
global J

nMeas = size(Z, 2);
% 实际上从这里的源代码可以看出
% 1st part: meas associated with existed target
mat1 = zeros(nMeas, maxTargIdx);
for i = 1 : nMeas                                   % 对每一个量测做假设，源于目标，新目标，杂波
    for j = 1 : maxTargIdx
        % find the targ cell whose index is j
        thisTarg = targ{cellfun(@(v) v{1}, targ) == j}; % 对v的操作是取第一个元素，这个操作针对于targ
        % 实际上从这里的源代码可以看出，
        % matTargIdx是递增，当一个目标的生命周期为0后，对于它的Cost值就变为无穷大，但是这个ID是被保留了的，那么这个目标thisTarg应该也不会删除
        if thisTarg{2} == 0 % lifePoint == 0, a disappeared targ
            mat1(i, j) = Inf;
        else
%             innov = H*thisTarg{3} - Z(:, i); % meas innovation 对每一个量测，取遍目标，做假设，然后量测循环
%             S = H*thisTarg{4} * H' + R; % prior cov of innov
%             dim = length(thisTarg{3});
%             x1 = log((1-pd)/pd); 
%             x2 = 0.5*(innov'*inv(S)*innov + dim*log(2*pi) + log(det(S)));   % 评分
%             mat1(i, j) = x1 + x2;
           %% 这里thisTarg是{idx lifePoint {100个采样}}
            mat1(i, j) = GenAverageMat1(thisTarg, Z(:,i),pd);
            mat1(i, j) = round(mat1(i, j),3);
        end
    end
end

% 2nd part: meas associated with new target
mat2 = inf(nMeas, nMeas); 
for i = 1 : nMeas
    %% 由于原始的程序densNew = 0 ,故
    mat2(i, i) = -log(densClt/20);%-log(densClt/20); %5是为了调试单个目标的程序，应该是20
    mat2(i, i) = round(mat2(i, i),3);
end

% 3rd part: meas associated with clutter
mat3 = inf(nMeas, nMeas);
for i = 1 : nMeas
    mat3(i, i) = -log(densClt/10); %10
    mat3(i, i) = round(mat3(i, i),3);
end

% put them together
ProbMat = [mat1, mat2, mat3];