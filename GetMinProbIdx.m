function chooseIdx = GetMinProbIdx(cellHypo,M,N,t)
% use the latest probabilities as the standard to judge
head = 2;
for i = 1 : N-1
    head = head + M^i;
end
rear = head + M^N - 1;      % ��һ�����У�Ԫ��22~85
arrayProb = cellfun(@(v) v{2}, cellHypo(head : rear));  % ��ȡprob
x = find(arrayProb == min(arrayProb)); % -log(prob), so find minimum ����С��-log(prob)Խ��
chooseBranch = ceil(x/M^(N-1));
if t == 3
    chooseBranch = 2;
end
idx = 1;
chooseIdx = [];

for i = 1 : N
    %% ����ֻȡ���һ�㱻ѡ�еļ�������
    if i == N           % ���һ��
        chooseIdx = [chooseIdx, (chooseBranch-1)*M^(i-1)+1 : chooseBranch*M^(i-1)];
%         chooseIdx = [chooseIdx, ...
%             idx+(chooseBranch-1)*M^(i-1)+1 : idx+chooseBranch*M^(i-1)];
       %% �ڱ����У���ΪM=4����������ͱ��64����������Χ��1��64
%         chooseIdx = chooseIdx - head;
    end
    idx = idx + M^i;
end
end