function zetaVarphi = getVarphi(Varphi,i)
%GETVARPHI 此处显示有关此函数的摘要
%   此处显示详细说明
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

