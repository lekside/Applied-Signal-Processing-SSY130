%% Question 1

clear all;
close all;
clc;

lc = lines(8);

ns = 16;
ne = 32;
ws = 25;

% NOISE INPUT
N = @(w) 0.1*(abs(w)>16).*(abs(w)<32);
w = linspace(-40,40,500);

% ANALOGIC FILTER
H = 1./(1+1i.*w./16);       % plot(w,abs(H)); grid on;
NH = N(w).*H;


% SAMPLING
wd = linspace(0,ws,1000);
k = -10:10;
NHd=zeros(size(k,2),size(wd,2));
for i=1:numel(k)
    NHd(i,:) = (ws/2/pi)*interp1(w,NH,wd+k(i)*ws,'linear','extrap');
end



%%%%%%%%% plots

X_index = find(any(sum(NHd,2),2));
k_index = X_index - max(k)-1;


helplines = [k_index; k_index(end)+1]*ws;
figure('Color','white','Position',[212.2000e+000   175.4000e+000   641.6000e+000   446.4000e+000]);
subplot(3,1,1)
plot(w,abs(NH), 'Linewidth',2); hold on;
plot( [helplines,helplines]', [zeros(size(helplines)),100*ones(size(helplines))]', 'k' )
ylim([0 1.1*max(max(abs(NH)))])
xlabel('w [\pi 10^{3} rad/s]')
ylabel('N(w)*H(w)')
title('N(w)*H(w)')


idx_kp = X_index( X_index > numel(k)/2 );
idx_kn = X_index( X_index < numel(k)/2 );
subplot(3,1,2)
plot(wd, abs(NHd(idx_kp,:)), '-','Color',lc(2,:),'Linewidth',2); hold on;
plot(wd, abs(NHd(idx_kn,:)), '-','Color',lc(4,:),'Linewidth',2);
plot([ws ws],[0 5*max(get(gca,'YLim'))], 'k' );
ylim([0 1.1*max(max(abs(NHd)))])
xlabel('w [\pi 10^{3} rad/s]')
ylabel('Nd(w)')
title('Nd(w) component-wise')
lgd=legend(num2str(k_index));
title(lgd,'k');


subplot(3,1,3)
plot(wd, sum(abs(NHd(X_index,:))),'Linewidth',2); hold on;
plot(wd, (ws/2/pi)*1/20.*wd./wd, 'Linewidth',1)
plot([ws ws],[0 5*max(max(abs(NHd)))], 'k' );
ylim([0 1.1*max(sum(abs(NHd(X_index,:)))) ])
xlabel('w [\pi 10^{3} rad/s]')
ylabel('Nd(w)')
title('Nd(w)')
legend({'Nd(w)','Threshold'})


set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'question1b-02'),'epsc')


%% Question 2

clear all;
close all;
clc;

f0 = 3;
fs = 10;

w0 = 2*pi*f0;
dt = 1/fs;
ws = 2*pi*fs;

w_dirac1 = @(k) +w0-k*ws;
w_dirac2 = @(k) -w0-k*ws;
Y_mag = @(w) abs( ws/(2*i) * exp(-i*pi*w/ws) * dt * sinc(w/ws) );

k=-3:3;

w= [];
mag = [];
for i=1:numel(k)
    w1 = w_dirac1(k(i));
    w2 = w_dirac2(k(i));
    w   = [w ; w1; w2];
    mag = [mag ; Y_mag(w1) ; Y_mag(w2) ];
end

figure('Color','white')
stem(w,mag,'filled');
hold on;
for i=1:numel(k)
    plot([k(i)*ws k(i)*ws], [0 max(get(gca,'YLim'))], '--','Color',0.7*[1 1 1]);
end
xlabel('w [10^{3} rad/s]')
ylabel('Y(w)')
title('Y(w) = H_{ZOH} * X_d(w)')
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'question2-01'),'epsc')


% answers
figure('Color','white')
ks = [-3 -2 -1 0 1 2 3];
stem(3-ks.*10, abs(sinc((3-ks.*10)/10)) , 'filled');
grid on;
xlabel('Sinus Frequency [kHz]')
ylabel('Sinus Magnitude')
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'question2-02'),'epsc')




