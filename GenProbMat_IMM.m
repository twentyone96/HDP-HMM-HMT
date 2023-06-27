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
      
function ProbMat = GenProbMat_IMM(maxTargIdx, targ, ...
    Z, H, R, pd, densNew, densClt)
global J
global IMMNUM
global dim

nMeas = size(Z, 2);

% 1st part: meas associated with existed target
mat1 = zeros(nMeas, maxTargIdx);
for i = 1 : nMeas                                   % 对每一个量测做假设，源于目标，新目标，杂波
    for j = 1 : maxTargIdx
        % find the targ cell whose index is j
        thisTarg = targ{cellfun(@(v) v{1}, targ) == j}; % 对v的操作是取第一个元素，这个操作针对于targ

%这里要要求mode probabilities 'im3ht algorthm a joint formulation of IMM and MHT'              
        cellIMM = thisTarg{3};
%        X = zeros(3*dim*IMMNUM, 1);
%        P = zeros(dim*3*IMMNUM,dim*3);
%        modeP = zeros(1,IMMNUM);
%        coeffi = zeros(1,IMMNUM);
%        zinv = zeros(dim,IMMNUM);
%        zcov = zeros(dim,IMMNUM);
        
	       X = cellIMM{1};
	       P = cellIMM{2};
	       coeffi = cellIMM{3};
	       modeP = cellIMM{4};

	       likelihood = zeros(IMMNUM,1);
	       for k=1:IMMNUM
	         innov = H*X(3*dim*(k-1)+1:3*dim*k) - Z(1:2, i); 
	         S = H*P(dim*3*(k-1)+1:dim*3*k,1:dim*3) * H.' + R; % prior cov of innov
             likelihood(k) = 0.5*(innov'*inv(S)*innov + dim*log(2*pi) + log(det(S)));   % 评分
             likelihood(k) = 1/likelihood(k);
%	         likelihood(k) = 1/(2*pi)+1/det(S)+exp(-0.5*(innov'*inv(S)*innov));   % 评分
	       end
	       modePPre = 0;
         for k=1:IMMNUM
           modePPre = modePPre + likelihood(k)*modeP(k);
	       end
	       for k=1:IMMNUM
           modeP(k) = likelihood(k)*modeP(k)/modePPre;
	       end
	       %%
	       coeffi = modeP;
        %%
         XX = zeros(3*dim,1);
         PP = zeros(dim*3,dim*3);
         for k=1:IMMNUM
             XX = XX + X(3*dim*(k-1)+1:3*dim*k)*modeP(k);
             PP = PP + P(dim*3*(k-1)+1:dim*3*k,1:dim*3)*modeP(k);
         end
         innov = H*XX - Z(1:2, i); % meas innovation 对每一个量测，取遍目标，做假设，然后量测循环
         S = H*PP * H' + R; % prior cov of innov
        if thisTarg{2} == 0 || (norm(Z(1:2, i)-H*XX,2) > 1000)% lifePoint == 0, a disappeared targ
            mat1(i, j) = Inf;
        else
             x1 = log((1-pd)/pd); 
             x2 = 0.5*(innov'*inv(S)*innov + dim*log(2*pi) + log(det(S)));   % 评分
             mat1(i, j) = x1 + x2;
           %% 这里thisTarg是{idx lifePoint {100个采样}}
%            mat1(i, j) = GenAverageMat1(thisTarg, Z(:,i),pd);
						cellIMM{1} = X;
						cellIMM{2} = P;
						cellIMM{3} = coeffi;
						cellIMM{4} = modeP;

		        thisTarg{3} = cellIMM;
        
        end
    end
end

% 2nd part: meas associated with new target
mat2 = inf(nMeas, nMeas); 
for i = 1 : nMeas
    %% 由于原始的程序densNew = 0 ,故
    mat2(i, i) = -log(densClt/20);%-log(densClt/20); %5是为了调试单个目标的程序，应该是20
end

% 3rd part: meas associated with clutter
mat3 = inf(nMeas, nMeas);
for i = 1 : nMeas
    mat3(i, i) = -log(densClt/10);
end

% put them together
ProbMat = [mat1/100, mat2, mat3];
