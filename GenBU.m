function [BU] = GenBU(T)
BU = zeros(9,3);
BU = [T^2/2  0   0;
       T     0   0;
       0     0   0;
       0   T^2/2 0;
       0   T     0;
       0   0     0;
       0   0     0;
       0   0     0;
       0   0     0];
end

