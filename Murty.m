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
[assgmt, cost] = Hungarian(costMat);    % 给出第一次mat用匈牙利算法求出的最优分配。
% 当漏检唯一的目标时，所有杂波对应的权重反而都是最小的，即所有量测点对应的都是杂波

assgmtTmp = assgmt;
assgmtTmp = [assgmtTmp;zeros(maxMes - size(assgmtTmp,1),1)];
solution{1} = {assgmtTmp, cost};           % 将分配保存

% xxxRec stands for 'record'
nodeRec = cell(1, 2);
assgmtRec = assgmt;

nodeList = MurtyPartition(nodeRec, assgmtRec, 1);   % 根据最优解求出murty算法需要的列表

while t < k                 % 求出后面所需的k-1个最优解
    tmp = []; % structure space for temporary (assgmt, cost) storage
    minCost = Inf;
    idxRec = -1;

    % try to find one node in the nodeList with the minimum cost 
    for i = 1 : size(nodeList, 2)
        node = nodeList{i};         % 分别提取列表中的每一个node
        Inclu = node{1};            % 提取include元素
        Exclu = node{2};            % 提取exclude元素
        mat = costMat;
        for j = 1 : size(Inclu, 1) % restrict: assignments must be included
            best = mat(Inclu(j, 1), Inclu(j, 2));   % 将最优解保留
            mat(Inclu(j, 1), :) = Inf;              % 删除最优解的一整行
            mat(Inclu(j, 1), Inclu(j, 2)) = best;   % 将最优解放入
        end
        for j = 1 : size(Exclu, 1) % restrict: assignments must be excluded
            mat(Exclu(j, 1), Exclu(j, 2)) = Inf;    % 将exlude的元素去除
        end
        
        [assgmt, cost] = Hungarian(mat);            % 求解处理之后的mat矩阵，求出此时的最优解
        
        if ismember(0, assgmt) % cost have to be inf
            continue;
        elseif cost < minCost                       % 代价最小的就是最优解
            minCost = cost;
            nodeRec = node;                         % node的定义见第40行，包含了include和exclude
            assgmtRec = assgmt;                     % 此时的最优解的对应关系
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
        solution{t} = {assgmTmp, minCost};         % 将上一个最优解存放
        idx = setdiff(1:size(nodeList, 2), idxRec);
        nodeList = [nodeList(idx), MurtyPartition(nodeRec, assgmtRec, 1)];  % 第一个参数node见第59行，node带入的include和exclude将影响list的元素
    end
end

