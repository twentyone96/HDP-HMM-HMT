function epsilonVarphi = getXi(Xi,i,ii)
%GETXI �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
global TL
global TargetNum
epsilonVarphi = zeros(2,1);
for targetIndex = 1 : TargetNum
    for j = 1 : TL
        for jj = 1 : TL
            if i == j && ii == jj
                epsilonVarphi(1) = epsilonVarphi(1) + Xi(i,ii,targetIndex);
            end
            if i == j && jj > ii
                epsilonVarphi(2) = epsilonVarphi(2) + Xi(i,ii,targetIndex);
            end
        end
    end
end
end

