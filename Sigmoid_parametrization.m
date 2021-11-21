function sigmoid_params = Sigmoid_parametrization(levels,g,ses)

f_sigmoid = @(x, g) (1-g) ./ (1 + exp(-x)) + g;
sigmoid_function = @(x,xdata)1+g./(1+exp(-(x(1)*xdata+x(2))));

total_sims = levels^2;
pars1 = [-3.5298; -3.6581; 5.0670]; % sigmoid model to go from synchronization to probability of correct response session1
pars2 = [-3.4476; -4.3652; 5.9553]; % sigmoid model to go from synchronization to probability of correct response session2
pars = [pars1 pars2];

scaling = linspace(1,1.5,levels);   % distance scaling (matches experimental design)
range = linspace(0.01,1,levels);    % contrast range (matches experimental design)
SCL = scaling' * ones(1, levels);
RNG = ones(levels, 1) * range;
v = [SCL(:),RNG(:), ones(total_sims, 1)];
% sigmoid_params = zeros(2,2);
emp_data_fit = f_sigmoid(v * pars(:,ses), g); % fitting BAT for 'ses'.
data = reshape(emp_data_fit,[20,20]);
x0 = [1,-1];
load(['k_BAT_ses',num2str(ses),'.mat']);
sigmoid_params = lsqcurvefit(sigmoid_function,x0,BAT,data);
