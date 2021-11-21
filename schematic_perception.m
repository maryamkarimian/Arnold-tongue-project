clear all
colormap jet
%%% settings %%%
N = 20;                             % oscillators placed along one dimension of the grid
Nn = N^2;                           % total number of oscillators
levels = 20;

dt = 1e-2;                          % integration time step
t_sim = 1;                          % simulation time (matches stimulus presentation time)
t_steps = floor(t_sim / dt) + 1;
t_average = floor(0.5 / dt):t_steps;

scaling = linspace(1,1.5,levels);   % distance scaling (matches experimental design)
range = linspace(0.01,1,levels);    % contrast range (matches experimental design)

%%% parameters %%%
diameter = 0.7;                     % annulus diameter (degree)
side_length = 7;                    % side length (degree) of square stimulus region
eccentricity = 7;                   % eccentricity (degree) of square stimulus region

kappa = 24.6301;                    % maximum coupling strength
lambda = 0.2227;                    % decay rate of coupling strength
offset = sqrt(eccentricity^2/2)-... % coordinate offset
    side_length/2; 

tic;
%%% functions %%%
F = @(theta, omega, N, C)...
    omega + 1 / N * sum(C .* sin(theta' - theta),2);

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

total_sims = levels^2;
BAT0 = zeros(levels);                
BAT1 = zeros(levels);                
BAT2 = zeros(levels);                
BAT3 = zeros(levels);


for r=1:10
    disp(r)
    bat0 = zeros(levels);
    bat1 = zeros(levels);
    bat2 = zeros(levels);
    bat3 = zeros(levels);
    coher = zeros(levels);
    for s=0:total_sims-1
        i = floor(s / levels) + 1;
        j = mod(s, levels) + 1;
        
        theta = zeros(Nn, t_steps);
        out_theta = zeros(Nn, t_steps);
        theta(:,1) = rand(Nn,1) * pi;
        out_theta(:,1) = rand(Nn,1) * pi;
        stimulus = generate_stimulus(scaling(i), range(end));
        outStim  = stimulus;
        stimulus = generate_stimulus(scaling(i), range(j));
        mean_luminance = mean(stimulus);
        mean_outLumin  = mean(outStim);
        norm_stim = (stimulus - mean_luminance).^2 / mean_luminance^2;
        norm_outStim = (outStim - mean_outLumin).^2 / mean_outLumin^2;
        
        contrast = sqrt(W * norm_stim);
        outCont  = sqrt(W * norm_outStim);
        freq = 25 * contrast + 25;
        outFreq = 25 * outCont + 25;
        omega = 2 * pi * freq;
        outOmega = 2 * pi * outFreq;
        for t=2:t_steps
            theta(:,t) = theta(:,t-1) + dt * F(theta(:,t-1), omega, Nn, K);
            out_theta(:,t) = out_theta(:,t-1) + dt * F(out_theta(:,t-1), outOmega, Nn, K);
        end
%         for t = t_average
%             theta_mean = mean(theta(:,t_average);
%             out_theta_mean = mean(out_theta(:,t_average));
%             coher(i,j) = coher(i,j) + abs(mean(exp([theta_mean ; out_theta_mean]*1i))); 
%         end
%         coher(i,j) = coher(i,j)/length(t_average);
        theta_mean = mean(mean(theta(:,t_average)));
        out_theta_mean = mean(mean(out_theta(:,t_average)));
        bat0(i,j) = mean(abs(mean(exp(theta(:,t_average)*1i))));
        bat1(i,j) = mean(abs(mean(exp(out_theta(:,t_average)*1i))));
        coher(i,j) = abs(mean(exp([theta_mean ; out_theta_mean]*1i)));
        bat3(i,j) = mean(abs(mean(exp(theta(:,t_average)*1i))))-coher(i,j);
%         bat4(i,j) = (bat0(i,j) - bat1(i,j)) / (bat0(i,j) + bat1(i,j));

    end 
    BAT0 = BAT0 + (bat0 - BAT0) / r;    % behavioral Arnold tongue (in_synch)
    BAT1 = BAT1 + (bat1 - BAT1) / r;    % behavioral Arnold tongue (out_synch)    
    BAT2 = BAT2 + (coher - BAT2) / r;    % BAT (synch(mean_thetha-mean_out_theta))
    BAT3 = BAT3 + (bat3 - BAT3) / r;    % BAT (in_synch - synch(mean_thetha-mean_out_theta))
%     imagesc(BAT,[0,1]);drawnow
end
toc;

dir1 = ['P:\FSE_MACSBIO\maryam.karimian\Arnold tongue project\'...
'experiments&simulations\simulations\learningSim_18\Discrimination\'];

figure;imagesc(BAT0)
xlabel('contrast'); ylabel('spacing');
p=colorbar;
colormap(lines(5)); col=colormap;
colormap(jet(255))
caxis([0 1])
% ylabel(p,'figure synchronization')
title('figure synchronization')
set(gcf,'color','w')
pbaspect([1 1 1])
nametif = [dir1,'FigureSynch.tif'];
namefig = [dir1,'FigureSynch.fig'];
saveas(gca, nametif);
saveas(gca, namefig);

figure;imagesc(BAT1)
xlabel('contrast'); ylabel('spacing');
p=colorbar;
colormap(lines(5)); col=colormap;
colormap(jet(255))
caxis([0 1])
% ylabel(p,'background synchronization')
title('ground synchronization')
set(gcf,'color','w')
pbaspect([1 1 1])
nametif = [dir1,'GroundSynch.tif'];
namefig = [dir1,'GroundSynch.fig'];
saveas(gca, nametif);
saveas(gca, namefig);

figure;imagesc(BAT2)
xlabel('contrast'); ylabel('spacing');
p=colorbar;
colormap(lines(5)); col=colormap;
colormap(jet(255))
caxis([0.4 0.8])
% ylabel(p,'background synchronization')
title('mean in-out phase coherence')
set(gcf,'color','w')
pbaspect([1 1 1])
nametif = [dir1,'in_out_synch.tif'];
namefig = [dir1,'in_out_synch.fig'];
saveas(gca, nametif);
saveas(gca, namefig);

figure;imagesc(BAT3)
xlabel('contrast'); ylabel('spacing');
p=colorbar;
colormap(lines(5)); col=colormap;
colormap(jet(255))
caxis([-1 1])
% ylabel(p,'background synchronization')
title('figure synchronization minus mean in-oit phase coherence')
set(gcf,'color','w')
pbaspect([1 1 1])
nametif = [dir1,'BATin_diff_InOutSynch.tif'];
namefig = [dir1,'BATin_diff_InOutSynch.fig'];
saveas(gca, nametif);
saveas(gca, namefig);


Discrimination.BAT_in = BAT0;
Discrimination.BAT_out = BAT1;
Discrimination.in_out_synch = BAT2;
Discrimination.BATin_diff_InOutSynch = BAT3;
fname = [dir1,'Discrimination.mat'];
save(fname, 'Discrimination');
