function stimulus = generate_stimulus(scaling,range)

omega = 5.7;                        % spatial frequency (cycles / degree)
diameter = 0.7;                     % annulus diameter (degree)
side_length = 7;                    % side length (degree) of square stimulus region
res = 50;                           % resolution of single grating
contrast_res = 480;                 % resolution of contrast image

%%% generate grating %%%
r = linspace(-diameter / 2,...
    diameter / 2,res);
[X,Y] = meshgrid(r);
radius = abs(X + Y*1i);
mask = radius <= diameter / 2;
grating = cos(2 * pi * radius * omega + pi);
grating = 0.5 * (grating.*mask + 1);

%%% generate stimulus %%%
new_diameter = scaling * diameter;
reps = ceil(side_length / new_diameter);
new_res = ceil(res * scaling);

stim_res = reps * new_res;
stimulus = ones(stim_res) * 0.5;

total = res^2;
for t=0:total-1
    i = floor(t / res);
    j = mod(t, res);
    contrast = rand() * range + 0.5 - range / 2;
    lower = 0.5 - contrast / 2;                % to keep the mean equals to 0.5.
    element = grating * contrast + lower;
    stimulus(i*new_res+1:i*new_res+res,...
        j*new_res+1:j*new_res+res) = element;
end

%%% remove outer border %%%

half_res = ceil((stim_res - contrast_res) / 2);
stimulus = stimulus(half_res+1:half_res+contrast_res,...
    half_res+1:half_res+contrast_res);

stimulus = stimulus(:);