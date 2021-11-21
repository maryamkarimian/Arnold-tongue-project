function [X,Y] = VF2Cort(X,Y)
% transform cortical coordinates to visual field coordinates using the
% inverse of the complex-logarithm transformation described in Schwartz
% (1980) - doi:10.1016/0042-6989(80)90090-5
%
% The function takes as input coordinates in the visual field and returns 
% coordinate on the cortex (with X and Y axes parallel and perpendicular 
% to the horizontal meridian, respectively).
%
% Parameter values are taken from:
% J.R. Polimeni, O.P. Hinds, M. Balasubramanian, A. van der Kouwe, 
% L.L. Wald, A.M. Dale, B. Fischl, E.L. Schwartz
% The human V1–V2–V3 visuotopic map complex measured via fMRI at 3 and 7 Tesla

k = 15.0;
a = 0.7;
b = 80;
alpha = .9; 

ecc = abs(X+Y*1i);
pa = angle(X+Y*1i);

Z = ecc.*exp(1i*alpha.*pa);
W = k*log((Z+a)./(Z+b))-k*log(a/b);

X = real(W);
Y = imag(W);