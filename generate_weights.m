function W = generate_weights(N, mu_x, mu_y, offset)
slope = 0.1720;                     % slope of rf size - eccentricity relationship 
intercept = -0.25;                  % intercept of rf size - eccentricity relationship 
theta = 1;                          % smallest allowed rf size 
conversion_factor = 1 / 1.7;        % sqrt(log(2))/sqrt(2), factr to convert between diameter and fwhm (based on Gaussian beams)
side_length = 7;                    % side length (degree) of square stimulus region
contrast_res = 480;                 % resolution of contrast image


f_diameter = @(x) max(slope * x + intercept, theta);

total_pix = contrast_res^2;
eccentricity = abs(mu_x + mu_y * 1i);
rf_diameter = f_diameter(eccentricity);
% fwhm = conversion_factor * rf_diameter;   % full width at half maximum
% sigma = fwhm / 2.355;                     %because fwhm = 2*sqrt(log(2))*sigma
% or simply: 
sigma = rf_diameter/4;                      % conversion_factor*2*sqrt(log(2))=4
r = linspace(0,side_length,contrast_res);
[X,Y] = meshgrid(r);
Y = flipud(Y);
X = X(:) + offset;
Y = Y(:) + offset;

W = zeros(total_pix, N);

for i=1:N
    R = (X - mu_x(i)).^2 + (Y - mu_y(i)).^2;
    W(:,i) = exp(- R ./ (2*sigma(i)^2));
end
D = diag(1./sum(W,1)); % semi-normalization (eq 2 in https://www.sciencedirect.com/science/article/pii/S0042698905005559?via%3Dihub)
W = (W * D)';

