function zetaVarphi = getVarphi(Varphi,i)
%GETVARPHI �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
global TL
zetaVarphi = zeros(2,1);
for j = 1 : TL
    for jj = 1 : TL 
        for jjj = 1 : TL 
            if jjj == i
                zetaVarphi(1) = zetaVarphi(1) + Varphi(j,jj,jjj);
            end
            if jjj > i
                zetaVarphi(2) = zetaVarphi(2) + Varphi(j,jj,jjj);
            end
        end
    end
end
end

