% KF_MHT_PREDICT  Perform Kalman Filter prediction step in MHT algorithm
%
% Syntax:
%   [X,P] = KF_PREDICT(X,P,A,Q,B,U)
%
% In:
%   cellTargSeed - a cell array containing description of targets. 
%       Each cell contains information of a full collection of targets
%       derived according to the corresponding hypo. One target was
%       described with {idx, lifePoint, X, P}. 
%   A - Transition matrix of discrete model (optional, default identity)
%   Q - Process noise of discrete model     (optional, default zero)
%   B - Input effect matrix                 (optional, default identity)
%   U - Constant input                      (optional, default empty)
%
% Out:
%  cellTargSeed - with X, P of each target in each cell modified
%   
% Description:
%   Perform Kalman Filter prediction step. The model is
%
%     x[k] = A*x[k-1] + B*u[k-1] + q,  q ~ N(0,Q).
% 
%   The predicted state is distributed as follows:
%   
%     p(x[k] | x[k-1]) = N(x[k] | A*x[k-1] + B*u[k-1], Q[k-1])
%
%   The predicted mean x-[k] and covariance P-[k] are calculated
%   with the following equations:
%
%     m-[k] = A*x[k-1] + B*u[k-1]
%     P-[k] = A*P[k-1]*A' + Q.
%
%   If there is no input u present then the first equation reduces to
%     m-[k] = A*x[k-1]

%KF_MHT_Predict_HDP(cellTargSeed_HDP{j}, F_CV,F_CA, Q_CV,Q_CA); 
%cellTarg_IMM的定义
%{id,life,
% {x1(2),p1(4),coefficient1(1),zinv1(2),zcov1(4)},
% {x2(2),p2(4),coefficient2(1),zinv2(2),zcov2(4)},
% {x3(2),p3(4),coefficient3(1),zinv3(2),zcov3(4)},
% {x4(2),p2(4),coefficient3(1),zinv4(2),zcov4(4)},}
function cellTargSeed = KF_MHT_Predict_IMM(cellTargSeed,...
    F_CV, F_CA, Q_CV, Q_CA,T)
global PTM;
global IMMNUM;
global dim 
global Q_coeff_IMM
if nargin < 2
    A = [];
end
if nargin < 3
    Q = [];
end
if nargin < 4
    B = [];
end
if nargin < 5
    U = [];
end
F = zeros(dim*3*IMMNUM,dim*3);
CTTemp9 = GenFCT(-9*pi/180,T);
CTTemp6 = GenFCT(6*pi/180,T);
F = [F_CV(1:dim*3,1:dim*3);F_CA(1:dim*3,1:dim*3);...
			CTTemp9(1:dim*3,1:dim*3);CTTemp6(1:dim*3,1:dim*3)];
Q = zeros(dim*3*IMMNUM,dim*3);
QTTemp9 = GenQCT(-9*pi/180,T);
QTTemp6 = GenQCT(6*pi/180,T);
Q = Q_coeff_IMM*[Q_CV(1:dim*3,1:dim*3);Q_CA(1:dim*3,1:dim*3);...
			QTTemp9(1:dim*3,1:dim*3);QTTemp6(1:dim*3,1:dim*3)];

for i = 1 : size(cellTargSeed, 2)
    oneCell = cellTargSeed{i};
    for j = 1 : size(oneCell, 2)
        % if lifePoint == 0, then pass without processing
        if oneCell{j}{2} == 0
            continue;
        end
        cellIMM = oneCell{j}{3};
%        X = zeros(3*dim*IMMNUM,1);
%        P = zeros(dim*3*IMMNUM,dim*3);
%        coeffi = zeros(1,IMMNUM);
%        modeP = zeros(1,IMMNUM);
%        zinv = zeros(dim,IMMNUM);
%        zcov = zeros(dim,IMMNUM);
        
	        X = cellIMM{1};
	        P = cellIMM{2};
	        coeffi = cellIMM{3};
	        modeP = cellIMM{4};

					%第j个模型
%	        coeffiPre = (PTM(:,:).*coeffi)/sum(PTM(:,:).*coeffi,2);
	        coeffiPre = zeros(IMMNUM,1);
	        coeffiPreMaxtrix = zeros(IMMNUM,IMMNUM);
	        for k=1:IMMNUM
	        	for kk=1:IMMNUM
		        	coeffiPre(k) = coeffiPre(k) + PTM(kk,k)*coeffi(kk);
	        	end
	        end
	        %%
	        modeP = coeffiPre;
	        %%
	        for k=1:IMMNUM
	        	for kk=1:IMMNUM
	        		coeffiPreMaxtrix(kk,k) = (PTM(kk,k)*coeffi(kk))/coeffiPre(k);
	        	end
	        end
	        xPre = zeros(3*dim*IMMNUM,1);
	        pPre = zeros(dim*3*IMMNUM,dim*3);
	        for k=1:IMMNUM
	        	for kk=1:IMMNUM
	        		xPre(3*dim*(k-1)+1:3*dim*k) = xPre(3*dim*(k-1)+1:3*dim*k) +...
	        			X(3*dim*(kk-1)+1:3*dim*kk)*coeffiPreMaxtrix(kk,k);
	        	end
	        end
	        for k=1:IMMNUM
	        	for kk=1:IMMNUM
	        		pPre(dim*3*(k-1)+1:dim*3*k,1:dim*3) = pPre(dim*3*(k-1)+1:dim*3*k,1:dim*3) +...
	        		P(dim*3*(kk-1)+1:dim*3*kk,1:dim*3) + coeffiPreMaxtrix(kk,k)*...
	        			(X(3*dim*(kk-1)+1:3*dim*kk)-xPre(3*dim*(k-1)+1:3*dim*k))*...
	        			(X(3*dim*(kk-1)+1:3*dim*kk)-xPre(3*dim*(k-1)+1:3*dim*k)).';
	        	end
	        end
	        for k=1:IMMNUM
	        	X(3*dim*(k-1)+1:3*dim*k) = F(dim*3*(k-1)+1:dim*3*k,1:dim*3)*...
	        		xPre(3*dim*(k-1)+1:3*dim*k);
	        	P(dim*3*(k-1)+1:dim*3*k,1:dim*3) = F(dim*3*(k-1)+1:dim*3*k,1:dim*3)*...
	           pPre(dim*3*(k-1)+1:dim*3*k,1:dim*3)*F(dim*3*(k-1)+1:dim*3*k,1:dim*3).'+...
							Q(dim*3*(k-1)+1:dim*3*k,1:dim*3);     		
            end
%            P = round(P,2);
            cellIMM{1} = X;
            cellIMM{2} = P;
            cellIMM{3} = coeffi;
            cellIMM{4} = modeP;

	        oneCell{j}{3} = cellIMM;
%         X = oneCell{j}{3};
%         P = oneCell{j}{4};
	        % Apply defaults
%	        if isempty(A)
%	            A = eye(size(X,1));
%	        end
%	        if isempty(Q)
%	            Q = zeros(size(X,1));
%	        end
%	        if isempty(B) && ~isempty(U)
%	            B = eye(size(X,1),size(U,1));
%	        end

	        % Perform prediction
%	        if isempty(U)
%	            X = A * X;
%	            P = A * P * A' + Q;
%	        else
%	            X = A * X + B * U;
%	            P = A * P * A' + Q;
%	        end

%	        oneCell{j}{3} = X;
%	        oneCell{j}{4} = P;
    end
    cellTargSeed{i} = oneCell;
end

