function [A, B, phy, yfit] = omSinefit( yvec, ang)
%% sine fit model:
% yvec = A*cos( 2* (ang - phy)) + B
%% test module: yvec = 0.5*cos(2*pi/10*(0:9)-pi/6)+3
% yvec = [3.0499, 3.3328, 3.4886, 3.4577, 3.2520, 2.9501, 2.6672, 2.5114, 2.5423, 2.7480];
%% 
% B作为拟合参数，效果很好
ang = ang / 180.0 * pi;
ang = ang(:);
fun = @( x, xdata) x(1) * cos( 2* (xdata - x(2))) + x(3);
options = optimset( 'TolFun', 1e-12, 'Display', 'off');
ang0L = find( yvec==max(yvec));
ang0 = ang(ang0L(1));
% x = lsqcurvefit( fun, [max(yvec)-min(yvec); 0; mean(yvec)], ang, yvec, [0; -pi; 0], [Inf; 0; Inf], options);
x = lsqcurvefit( fun, [max(yvec)-min(yvec); ang0; mean(yvec)], ang, yvec, [0; 0; 0], [Inf; pi; Inf], options);
A = x(1); phy = x(2); B = x(3);
%% evaluate
yfit = A * cos( 2* (ang - phy)) + B;
% phy need to be in [0,pi)
phy = mod( phy, pi);
phy = phy / pi * 180;











