clear all
%%% settings %%%
N = 20;                             % oscillators placed along one dimension of the grid
Nn = N^2;                           % total number of oscillators
levels = 20;
total_sims = levels^2;
dt = 1e-2;                          % integration time step
t_sim = 1;                          % simulation time (matches stimulus presentation time)
t_steps = floor(t_sim / dt) + 1;
t_average = floor(0.5 / dt):t_steps;

scaling = linspace(1,1.5,levels);   % distance scaling (matches experimental design)
range = linspace(0.01,1,levels);    % contrast range (matches experimental design)

%%% parameters %%%
g = .5;
diameter = 0.7;                     % annulus diameter (degree)
side_length = 7;                    % side length (degree) of square stimulus region
eccentricity = 7;                   % eccentricity (degree) of square stimulus region
kappa = 24.6301;                    % maximum coupling strength
lambda = 0.2227;                    % decay rate of coupling strength
EPSILON = 3.87;%7.67;%3.069;              % learning rate
offset = sqrt(eccentricity^2/2)-side_length/2; % coordinate offset 

tic;
%%% functions %%%
F = @(theta, omega, N, C)...
    omega + 1 / N * sum(C .* sin(theta' - theta),2);

sigmoid_function = @(x,xdata)g+1./(1+exp(-(x(1)*xdata+x(2))));

%%% initializations %%%
r = linspace(0, side_length, N);
[X, Y] = meshgrid(r);
Y = flipud(Y);
X = X(:) + offset;                  % visual field x-coordinates
Y = Y(:) + offset;                  % visual field y-coordinates
W = generate_weights(Nn,X,Y,offset);
[Xc,Yc] = VF2Cort(X,Y);             % cortex x- and y-coordinates
CD = sqrt((Xc - Xc').^2 +...        % pairwise cortical distances
    (Yc - Yc').^2);
K = kappa * exp(-lambda * CD);      % pairwise coupling strength 
sigmoid_params = Sigmoid_parametrization(levels,g,1);

%% Session 1

ses = 1;
disp('ses1');
%%% initializations %%%
BAT = zeros(levels);                % behavioral Arnold tongue
Q = zeros(Nn);                      % "correlation" / FC matrix
cosmat = zeros(Nn,Nn,t_steps);      % matrix of cos(theta_i - theta_j)
for r=1:3
    disp(r)
    bat = zeros(levels);
    q = zeros(Nn);
    for s=0:total_sims-1
        i = floor(s / levels) + 1;
        j = mod(s, levels) + 1;
        theta = zeros(Nn, t_steps);
        theta(:,1) = rand(Nn,1) * pi;
        stimulus = generate_stimulus(scaling(i), range(j));
        mean_luminance = mean(stimulus);
        norm_stim = (stimulus - mean_luminance).^2 / mean_luminance^2;
        contrast = sqrt(W * norm_stim);
        freq = 25 * contrast + 25;
        omega = 2 * pi * freq;
        for t=2:t_steps
            theta(:,t) = theta(:,t-1) + dt * F(theta(:,t-1), omega, Nn, K);
            cosmat(:,:,t) = cos(theta(:,t-1) - theta(:,t-1)');
        end
        bat(i,j) = mean(abs(mean(exp(theta(:,t_average)*1i))));
        prob = sigmoid_function(sigmoid_params,bat(i,j)); 
        mask = cosmat >=0;
        cosmat = cosmat .* mask; 
        q = q + prob * (mean(cosmat(:,:,t_average),3) - q) / (s+1);
    end 
    BAT = BAT + (bat - BAT) / r;    % behavioral Arnold tongue (in_synch)
    Q = Q + (q - Q) / r;
end
fname = ['eps1_k_AT_ses',num2str(ses),'.mat'];
save(fname,'K','BAT','Q');

%% session 2-8

for ses = 2:8
    disp('ses')
    disp (ses)
    %%% initializations %%%
    K = K.*exp(-EPSILON) + (1-exp(-EPSILON)) * kappa * Q;
    BAT = zeros(levels);                % behavioral Arnold tongue
    Q = zeros(Nn);                      % "correlation" / FC matrix
    cosmat = zeros(Nn,Nn,t_steps);      % matrix of cos(theta_i - theta_j)
    
    for r=1:10
        disp(r)
        bat = zeros(levels);
        q = zeros(Nn);
        for s=0:total_sims-1
            i = floor(s / levels) + 1;
            j = mod(s, levels) + 1;
            theta = zeros(Nn, t_steps);
            theta(:,1) = rand(Nn,1) * pi;
            stimulus = generate_stimulus(scaling(i), range(j));
            mean_luminance = mean(stimulus);
            norm_stim = (stimulus - mean_luminance).^2 / mean_luminance^2;
            contrast = sqrt(W * norm_stim);
            freq = 25 * contrast + 25;
            omega = 2 * pi * freq;
            for t=2:t_steps
                theta(:,t) = theta(:,t-1) + dt * F(theta(:,t-1), omega, Nn, K);
                cosmat(:,:,t) = cos(theta(:,t-1) - theta(:,t-1)');
            end
            bat(i,j) = mean(abs(mean(exp(theta(:,t_average)*1i))));
            prob = sigmoid_function(sigmoid_params,bat(i,j));
            mask = cosmat >=0;
            cosmat = cosmat .* mask;
            q = q + prob * (mean(cosmat(:,:,t_average),3) - q) / (s+1);
        end
        BAT = BAT + (bat - BAT) / r;    % behavioral Arnold tongue (in_synch)
        Q = Q + (q - Q) / r;
    end
    fname = ['eps1_k_AT_ses',num2str(ses),'.mat'];
    save(fname,'K','BAT','Q');
end