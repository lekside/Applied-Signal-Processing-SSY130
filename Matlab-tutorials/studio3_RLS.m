%% 
% Demo of RLS filtering
% 
clear all
%close all
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
x = x0 +  0.5* randn(N,1);
% This is the desired signal

% Create space to save filter coefficients
hhat = zeros(Nh,N); % estimated state 
e = zeros(N,1);

Rinv = 10*eye(10); %Initalization
alpha = 1;
% step size
for k=Nh:N,
  e(k) = x(k) - hhat(:,k)'*y(k-Nh+1:k);
  K = Rinv*y(k-Nh+1:k);
  Rinv = (Rinv - K*K'/(alpha + y(k-Nh+1:k)'*K))/alpha;
  hhat(:,k+1) = hhat(:,k) + Rinv*(y(k-Nh+1:k)*e(k));
end

% Estimate residual variance
% calculated based on last 100 samples
sigma_e2 = sum(e((end-99):end).^2)/100


% Plot filter error over time
figure(1);
pp=plot(e);
set(pp,'LineWidth',2)
p = gca;
set(p,'FontSize',14)
title(['e(n) RLS \alpha=',num2str(alpha),'  Var(e(n)) =', num2str(sigma_e2)]);

%  Plot estimated filter coefficients over time

figure(2)
pp=plot(hhat');
set(pp,'LineWidth',2);
p = gca;
set(p,'FontSize',14);
title(['hhat(n) RLS \alpha=',num2str(alpha)]);

%% Now add forgetting to the algorihtm and see the effect

hhat = zeros(Nh,N); % clear estimated filter values 

Rinv = 100*eye(10); %Initalization
% Forgetting factor
alpha = 0.99;
for k=Nh:N
  e(k) = x(k) - hhat(:,k)'*y(k-Nh+1:k);
  K = Rinv*y(k-Nh+1:k);
  Rinv = (Rinv - K*K'/(alpha + y(k-Nh+1:k)'*K))/alpha;
  hhat(:,k+1) = hhat(:,k) + Rinv*y(k-Nh+1:k)*e(k);
end

% Estimate residual variance
% calculated based on last 100 samples
sigma_e2 = sum(e(end-Nh-100:end-Nh-1).^2)/100


% Plot filter error over time
figure(3);
pp=plot(e);
set(pp,'LineWidth',2)
p = gca;
set(p,'FontSize',14)
title(['e(n) RLS \alpha=',num2str(alpha),'  Var(e(n)) =', num2str(sigma_e2)]);


%  Plot estimated filter coefficients over time
figure(4)
pp=plot(hhat');
set(pp,'LineWidth',2);
p = gca;
set(p,'FontSize',14)
title(['hhat(n) RLS \alpha=',num2str(alpha)]);



%% Things to test and think about
% Test alpha 0.8, 0.9 0.99 
% How is the convergence speed is changed?
% How is the residual variance changed?