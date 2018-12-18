%% Demonstration of Kalman filtering
% Define the motion model and generate a movement
T = 0.1;
tvec = 0.1*(0:2999);
A = [ 1 T ; 0 1]; % first state is position, second state is velocity
B = [0 1]';
U = [5; zeros(999,1); -5; zeros(999,1); 5; zeros(999,1)]+ 0.00*randn(3000,1);
X = ltitr(A,B,U); %
subplot(211);
plot(tvec,[X(:,1) ], '-','linewidth', 2); ylabel('position [m]','fontsize',16); xlabel('time [s]')
p = gca();
set(p,'fontsize',16)
subplot(212);
plot(tvec,[X(:,2) ], '-','linewidth', 2); ylabel('velocity [m]'); xlabel('time [s]')
p = gca();
set(p,'fontsize',16)

%% Measurement
C = [1 0];
y = X*C'+1*randn(3000,1);
subplot(111);
plot(tvec,[y(:,1) ], '-','linewidth', 2); ylabel('position [m]','fontsize',16); xlabel('time [s]')
title('position measurement')
p = gca();
set(p,'fontsize',16)

%% See the trivial speed estimator
plot(tvec(1:end-1),[  diff(y)/T X(1:end-1,2)],'linewidth', 2); ylabel('velocity [m]'); xlabel('time [s]')
hand = legend('Diff estimate', 'True speed');
%set(hand,'fontsize',14);
p = gca();
set(p,'fontsize',16)
%% Form the Kalman filter

Q = [0  0; 0 1]; % State-noise covariance
R = 1e3;         % Measurment noise covariance (try 1e-3, 1, 1e3 and 1e6) 

P = 1e3*eye(2);   % Initial state-covaraince
hx = zeros(2,3000); % estimated state 
hxp = zeros(2,3000); % measurement updated state


for k=1:3000-1,
  hxp(:,k) = hx(:,k) + P*C'*inv(C*P*C'+R)*(y(k)-C*hx(:,k)); % measurement update
  Pp = P - P*C'*inv(C*P*C'+R)*C*P; 
  hx(:,k+1) = A * hxp(:,k);  % State-update (time update)
  P = A*Pp*A' + Q;
end
hxp(:,3000) = hx(:,3000) + P*C'*inv(C*P*C'+R)*(y(3000)-C*hx(:,3000));


%% Position estimates
plot([  X(:,1)-y, X(:,1)-hxp(1,:).' ],'linewidth', 2)
%plot([X(:,1)-hxp(1,:)' ])
hand = legend('Measurement error','Kalman estimate error' );
set(hand,'fontsize',14);
p = gca();
set(p,'fontsize',16);


%% Position estimates
% plot(cumsum(abs([X(:,1)-hxp(1,:)' X(:,1)-y])),'linewidth', 2)
% hand = legend('Cummulative Kalman abs error', 'Cummulative measurement abs error')
% set(hand,'fontsize',14);


%% Speed estimates
figure(1);
plot([hxp(2,1:end)', X(:,2)  ],'linewidth', 2)
hand= legend('Kalman estimated speed', 'True speed' ,'Location','southeast');
%set(hand,'fontsize',14);
title(['R=1e', num2str(log10(R))],'fontsize', 14);
%print -dpng R5.png
p = gca();
set(p,'fontsize',16);




