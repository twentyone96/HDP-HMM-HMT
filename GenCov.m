% Bayesian approach extended object and cluster tracking using Random
% Matrices, apex E provides the related procedures
% translate cov in POL to cov in Cart
function [a2_R] = GenCov(x_stateC)
%GENMEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here
  global RadarLocX
  global RadarLocY
  global CovA
  global re_d
  global CartX
  global CartY
  % cov of x_state
%  a2_R= zeros(2,2);
  a2_R= zeros(3,3);
  CovA2=(CovA)^2;
  CovB2=(re_d)^2;
  %这里的参数要与GenMeasurmentsPol统一
%   % all pairs in list has been removed
% 

    [theta,rho] = cart2pol(x_stateC(1)-RadarLocX,x_stateC(2)-RadarLocY);
    theta = theta + pi/2;
    D2 = rho^2;
     a2_R(1,1) = CovB2*(sind(theta))^2+D2*(cosd(theta))^2*CovA2;
     a2_R(1,2) = (CovB2 - D2*CovA2)*cosd(theta)*sind(theta);
     a2_R(2,1) = a2_R(1,2);
     a2_R(2,2) = CovB2*(cosd(theta))^2+D2*(sind(theta))^2*CovA2;

   a2_R(1,1) = CartX^2;
   a2_R(1,2) = 0;
   a2_R(2,1) = 0;
   a2_R(2,2) = CartY^2;
end