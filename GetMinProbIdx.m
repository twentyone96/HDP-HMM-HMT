function chooseIdx = GetMinProbIdx(cellHypo,M,N,t)
% use the latest probabilities as the standard to judge
head = 2;
for i = 1 : N-1
    head = head + M^i;
end
rear = head + M^N - 1;      % 第一次运行，元素22~85
arrayProb = cellfun(@(v) v{2}, cellHypo(head : rear));  % 提取prob
x = find(arrayProb == min(arrayProb)); % -log(prob), so find minimum 概率小，-log(prob)越大
chooseBranch = ceil(x/M^(N-1));
if t == 3
    chooseBranch = 2;
end
idx = 1;
chooseIdx = [];

for i = 1 : N
    %% 这里只取最后一层被选中的假设索引
    if i == N           % 最低一层
        chooseIdx = [chooseIdx, (chooseBranch-1)*M^(i-1)+1 : chooseBranch*M^(i-1)];
%         chooseIdx = [chooseIdx, ...
%             idx+(chooseBranch-1)*M^(i-1)+1 : idx+chooseBranch*M^(i-1)];
       %% 在本例中，因为M=4，到第三层就变成64，故索引范围从1：64
%         chooseIdx = chooseIdx - head;
    end
    idx = idx + M^i;
end
end