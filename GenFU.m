function [FU] = GenFU(T)
FU = zeros(9,9);
FU = [ 1 T T^2/2 0  0  0  0  0  0;
      0  1  T  0  0  0  0  0  0;
      0  0  1  0  0  0  0  0  0;
      0  0  0  1  T T^2/2 0  0  0;
      0  0  0  0  1  T  0  0  0;
      0  0  0  0  0  1  0  0  0;
      0  0  0  0  0  0  0  0  0;
      0  0  0  0  0  0  0  0  0;
      0  0  0  0  0  0  0  0  0];
end

