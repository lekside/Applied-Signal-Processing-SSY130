close all;
clc;
if ~ exist(fullfile(pwd,'images'),'dir'), mkdir images; end


%% Question 2

f_s = 16000;    % Hz
h = flipud(h_reversedOrder_1);
[H,w] = freqz(h,1,512,f_s); %'whole',


figure('Color','white')
stem(0:numel(h)-1, h); grid on;
xlabel 'Samples', ylabel 'h(n) - Filter coefficients', title 'Estimated Filter h(n)'
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'images/','h-coeffs'),'epsc')

figure('Color','white');
subplot(2,1,1);
plot(w/f_s,20*log10(abs(H))); grid on;
xlabel 'Normalized Frequency (\times \pi [rad/sample])', ylabel 'Magnitude (dB)', title 'Estimated Filter H(w)'
subplot(2,1,2)
plot(w/f_s, 180/pi*wrapToPi(angle(H))); grid on;
xlabel 'Normalized Frequency (\times \pi [rad/sample])', ylabel 'Phase (degrees)'
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'images/','H'),'epsc')





%% Question 3

f_s = 16000;    % Hz
h = [1 0 0 0.5];
[H,w] = freqz(h,1,512,f_s); %'whole',

figure('Color','white');
subplot(2,1,1);
plot(w/f_s,20*log10(abs(H))); grid on;
xlabel 'Normalized Frequency (\times \pi [rad/sample])', ylabel 'Magnitude (dB)', title 'Frequency response H(w) of h = [1, 0, 0, 0.5]'
subplot(2,1,2)
plot(w/f_s, 180/pi*wrapToPi(phase(H))); grid on;
xlabel 'Normalized Frequency (\times \pi [rad/sample])', ylabel 'Phase (degrees)'
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'images/','H-echo'),'epsc')


%% Question 6


f_s = 16000;    % Hz

h = flipud(h_reversedOrder_max);
figure('Color','white')
stem(0:numel(h)-1, h); grid on;
xlabel 'Samples', ylabel 'h(n) - Filter coefficients', title 'Estimated Filter h(n)'
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'images/','h-coeffs-SIN-MAX'),'epsc')

h = flipud(h_reversedOrder_100);
figure('Color','white')
stem(0:numel(h)-1, h); grid on;
xlabel 'Samples', ylabel 'h(n) - Filter coefficients', title 'Estimated Filter h(n)'
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'images/','h-coeffs-SIN-100'),'epsc')

h = flipud(h_reversedOrder_10);
figure('Color','white')
stem(0:numel(h)-1, h); grid on;
xlabel 'Samples', ylabel 'h(n) - Filter coefficients', title 'Estimated Filter h(n)'
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'images/','h-coeffs-SIN-10'),'epsc')

[H,w] = freqz(h,1,512,f_s); %'whole',
figure('Color','white');
subplot(2,1,1);
plot(w/f_s,20*log10(abs(H))); grid on;
xlabel 'Normalized Frequency (\times \pi [rad/sample])', ylabel 'Magnitude (dB)', title 'Estimated Filter H(w)'
subplot(2,1,2)
plot(w/f_s, 180/pi*wrapToPi(angle(H))); grid on;
xlabel 'Normalized Frequency (\times \pi [rad/sample])', ylabel 'Phase (degrees)'
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'images/','H-SIN-MAX'),'epsc')

%%%%%%%%%%%%

h = flipud(h_reversedOrder_MAX_inc);
figure('Color','white')
stem(0:numel(h)-1, h); grid on;
xlabel 'Samples', ylabel 'h(n) - Filter coefficients', title 'Estimated Filter h(n)'
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'images/','h-coeffs-SIN-MAX-inc'),'epsc')

h = flipud(h_reversedOrder_100_inc);
figure('Color','white')
stem(0:numel(h)-1, h); grid on;
xlabel 'Samples', ylabel 'h(n) - Filter coefficients', title 'Estimated Filter h(n)'
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'images/','h-coeffs-SIN-100-inc'),'epsc')

h = flipud(h_reversedOrder_10_inc);
figure('Color','white')
stem(0:numel(h)-1, h); grid on;
xlabel 'Samples', ylabel 'h(n) - Filter coefficients', title 'Estimated Filter h(n)'
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'images/','h-coeffs-SIN-10-inc'),'epsc')
