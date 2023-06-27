function [F_CT] = GenFCT(turnRate,T)
F_CT = zeros(9,9);
% F_CT(1,:) = [1, sin(turnRate*T)/turnRate, (1-cos(turnRate*T))/turnRate^2,0,0,0,0,0,0];
% F_CT(2,:) = [0, cos(turnRate*T), sin(turnRate*T)/turnRate,0,0,0,0,0,0];
% F_CT(3,:) = [0, -turnRate*sin(turnRate*T), cos(turnRate*T),0,0,0,0,0,0];
% F_CT(4,:) = [0,0,0,1, sin(turnRate*T)/turnRate, (1-cos(turnRate*T))/turnRate^2,0,0,0];
% F_CT(5,:) = [0,0,0,0, cos(turnRate*T), sin(turnRate*T)/turnRate,0,0,0];
% F_CT(6,:) = [0,0,0,0, -turnRate*sin(turnRate*T), cos(turnRate*T),0,0,0];
% F_CT(7,:) = [0,0,0,0,0,0,1, sin(turnRate*T)/turnRate, (1-cos(turnRate*T))/turnRate^2];
% F_CT(8,:) = [0,0,0,0,0,0,0, cos(turnRate*T), sin(turnRate*T)/turnRate];
% F_CT(9,:) = [0,0,0,0,0,0,0, -turnRate*sin(turnRate*T), cos(turnRate*T)];
F_CT=[1,(sin(turnRate*T))/turnRate,0,0,-(1-cos(turnRate*T))/turnRate,0,0,0,0;
   0,cos(turnRate*T),0,0,-sin(turnRate*T),0,0,0,0;
   0,0,1,0,0,0,0,0,0;
   0,(1-cos(turnRate*T))/turnRate,0,1,sin(turnRate*T)/turnRate,0,0,0,0;
   0,sin(turnRate*T),0,0,cos(turnRate*T),0,0,0,0;
   0,0,0,0,0,1,0,0,0;
   0,0,0,0,0,0,1,0,0;
   0,0,0,0,0,0,0,1,0;
   0,0,0,0,0,0,0,0,1];    % ÈıÎ¬µÑ¿¨¶û×ø±êÏµ
end

