% Murty's algorithm finds out the kth minimum assignments, k = 1, 2, ...
% 
% Syntax: 
%   solution = Murty(costMat, k)
%   
% In: 
%   costMat - nMeas*nTarg cost matrix.
%   k - the command number controlling the output size.
%   
% Out: 
%   solution - cell array containing the minimum, 2nd minimum, ..., 
%       kth minimum assignments and their costs. Each solution{i} 
%       contains {assgmt, cost} where assgmt is an nMeas*1 matrix 
%       giving the ith minimum assignment; cost is the cost of this 
%       assignment.
      
function solution = Murty(costMat, k)

global maxMes

% asso = zeros(maxMes,1);
% prob = zeros(1,1);
% hypo = {asso, prob}
solution = cell(1, k);
%solution{1:k} = {hypo};

t = 1;
[assgmt, cost] = Hungarian(costMat);    % ������һ��mat���������㷨��������ŷ��䡣
% ��©��Ψһ��Ŀ��ʱ�������Ӳ���Ӧ��Ȩ�ط���������С�ģ�������������Ӧ�Ķ����Ӳ�

assgmtTmp = assgmt;
assgmtTmp = [assgmtTmp;zeros(maxMes - size(assgmtTmp,1),1)];
solution{1} = {assgmtTmp, cost};           % �����䱣��

% xxxRec stands for 'record'
nodeRec = cell(1, 2);
assgmtRec = assgmt;

nodeList = MurtyPartition(nodeRec, assgmtRec, 1);   % �������Ž����murty�㷨��Ҫ���б�

while t < k                 % ������������k-1�����Ž�
    tmp = []; % structure space for temporary (assgmt, cost) storage
    minCost = Inf;
    idxRec = -1;

    % try to find one node in the nodeList with the minimum cost 
    for i = 1 : size(nodeList, 2)
        node = nodeList{i};         % �ֱ���ȡ�б��е�ÿһ��node
        Inclu = node{1};            % ��ȡincludeԪ��
        Exclu = node{2};            % ��ȡexcludeԪ��
        mat = costMat;
        for j = 1 : size(Inclu, 1) % restrict: assignments must be included
            best = mat(Inclu(j, 1), Inclu(j, 2));   % �����ŽⱣ��
            mat(Inclu(j, 1), :) = Inf;              % ɾ�����Ž��һ����
            mat(Inclu(j, 1), Inclu(j, 2)) = best;   % �����Ž����
        end
        for j = 1 : size(Exclu, 1) % restrict: assignments must be excluded
            mat(Exclu(j, 1), Exclu(j, 2)) = Inf;    % ��exlude��Ԫ��ȥ��
        end
        
        [assgmt, cost] = Hungarian(mat);            % ��⴦��֮���mat���������ʱ�����Ž�
        
        if ismember(0, assgmt) % cost have to be inf
            continue;
        elseif cost < minCost                       % ������С�ľ������Ž�
            minCost = cost;
            nodeRec = node;                         % node�Ķ������40�У�������include��exclude
            assgmtRec = assgmt;                     % ��ʱ�����Ž�Ķ�Ӧ��ϵ
            idxRec = i;
        end
    end
    
    
    if idxRec == -1 % all node in the nodeList leads to inf cost
        for i = t+1 : k
            solution{i} = solution{t};
        end
        t = k;
    else    
        t = t + 1;
        assgmTmp = assgmtRec;
		assgmTmp = [assgmTmp;zeros(maxMes - size(assgmTmp,1),1)];
        solution{t} = {assgmTmp, minCost};         % ����һ�����Ž���
        idx = setdiff(1:size(nodeList, 2), idxRec);
        nodeList = [nodeList(idx), MurtyPartition(nodeRec, assgmtRec, 1)];  % ��һ������node����59�У�node�����include��exclude��Ӱ��list��Ԫ��
    end
end

