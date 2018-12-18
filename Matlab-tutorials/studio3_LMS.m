%% 
% Demo of LMS filtering
%% 
N = 1600;
% Start with a FIR filter
h0 = (10:-1:1)';
%h1 = h0;
h1 = h0-2; % System to change to at time = 400
Nh = length(h0);

% Use white noise as input signal
y = randn(N,1);

x0 = zeros(N,1);
% create output signal
for k=Nh:N
  if k<400
    x0(k) = h0'*y(k-Nh+1:k);
  else % system change after k=400
    x0(k) = h1'*y(k-Nh+1:k);
  end
end

% Add measurement noise
x = x0 + 0.1 * randn(N,1);
% This is the desired signal

% Create space to save filter coefficients
hhat = zeros(Nh,N); % clear estimated filter 
e = zeros(N,1);

% step size
mu = 2e-2;
for k=Nh:(N)
  e(k) = x(k) - hhat(:,k)'*y((k-Nh+1):k);
  hhat(:,k+1) = hhat(:,k) + 2*mu*y((k-Nh+1):k)*e(k);
end
figure(1)

% Estimate residual variance
% calculated based on last 100 samples
sigma_e2 = sum(e((end-99):end).^2)/100

% Plot filter error over time
pp=plot(e);
set(pp,'LineWidth',2)
p = gca;
set(p,'FontSize',14)
title(['e(n) LMS \mu=',num2str(mu),'  Var(e(n)) =', num2str(sigma_e2)]);

%  Plot estimated filter coefficients over time
figure(2)
pp=plot(hhat')
set(pp,'LineWidth',2)
p = gca;
set(p,'FontSize',14)
title(['hhat(n) LMS \mu=',num2str(mu)]);

%% Things to test and think about
% 1) what  is the expected residualt variance (i.e. variance of e)
% when the filter has converged?
%
% 2) Try shorter and longer step lenghts. How is the convergence speed
% affected?
% 3) Is the resiudual variance affected by different values of the step
% length?

% Now try the same setup but with a different input.


%%
y = sin(0.05*pi*(1:N)') + 0*randn(N,1); 
% You can test with different levels of the randn component, i.e. 0, 0.1, 1 

x0 = zeros(N,1);
% create output signal
for k=Nh:N
  if k<400
    x0(k) = h0'*y(k-Nh+1:k);
  else % system change after k=400
    x0(k) = h1'*y(k-Nh+1:k);
  end
end

% Add measurement noise
% This is the desired signal 
x = x0 + 0.1 * randn(N,1);

% Create space to save filter coefficients
hhat = zeros(Nh,N); % clear variable
e = zeros(N,1);

% step size
mu = 1e-2;
for k=Nh:N
  e(k) = x(k) - hhat(:,k)'*y((k-Nh+1):k);
  hhat(:,k+1) = hhat(:,k) + mu*y((k-Nh+1):k)*e(k);
end
figure(1)

% Estimate residual variance
% calculated based on last 100 samples
sigma_e2 = sum(e((end-99):end).^2)/100

% Plot filter error over time
pp=plot(e);
set(pp,'LineWidth',2)
p = gca;
set(p,'FontSize',14)
title(['e(n) LMS \mu=',num2str(mu),'  Var(e(n)) =', num2str(sigma_e2)]);

%  Plot estimated filter coefficients over time
figure(2)
pp=plot(hhat')
set(pp,'LineWidth',2)
p = gca;
set(p,'FontSize',14)
title(['hhat(n) LMS \mu=',num2str(mu)]);

% Did the error e converge towards zero?

% Did the filter coefficients converge to their true values?
% Why not?