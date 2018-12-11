% Params
f_s = 16000;            % Hz
w_avg = 2*pi*(4000);    % rad/s
v_sound = 343;          % m/s
L=8;                    % up-sampling factor
char2symbols = 8/2;     % 1 char = 8 bits;  2 bits = 1 symbol
pilot_chars = 16+8;     % 16 pilot + 8 cp

N_samples = pilot_chars * char2symbols *L;  % number of transmitted samples

t0 = @(theta) theta/w_avg;
ds = @(theta) t0(theta) * v_sound;
dt = N_samples / f_s;

v = @(theta) ds(theta) / dt;    % m/s

% Moving to towards DSP (theta=20deg)
v(20/180*pi) *100   % cm/s
%     ans =
% 
%         9.9248

% Moving further away from DSP (theta=-15deg)
v(-15/180*pi) *100  % cm/s
%     ans =
% 
%        -7.4436



% Worst case - change w_avg to 5000 Hz
v(45/180*pi) *100  % cm/s
%     ans =
% 
%        17.8646