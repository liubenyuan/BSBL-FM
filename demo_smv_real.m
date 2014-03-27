% This demo shows the capability of recoverying REAL valued signals
% using the proposed algorithm : BSBL-FM
%
% author: liubenyuan@gmail.com
% date:   2013-03-04
%
clear all;  close all;
%==========================================================================
rng(1985,'v4');

% problem dimension
M      = 128;     % row number of the dictionary matrix 
N      = 256;     % column number
blkNum = 5;       % nonzero block number
blkLen = 16;      % block length
iterNum= 1;       % number of experiments (100)

% Generate the known matrix with columns draw uniformly from the surface of a unit hypersphere
Phi = randn(M,N);
for i=1:N
    Phi(:,i) = Phi(:,i) / norm(Phi(:,i));
end

% load data
load demo.mat;

% prepare
Wgen     = re; 
% compressed the signal
signal   = Phi * Wgen;
% Observation noise
SNR      = 15;
stdnoise = std(signal)*10^(-SNR/20);
noise    = randn(M,1) * stdnoise;
% Noisy observation
Y        = signal + noise;

%=== BSBL-FM ==============================================================
blkStartLoc = 1:blkLen:N;
learnLambda = 1;

tic;
Result2 = BSBL_FM(Phi,Y,blkStartLoc,learnLambda,'learnType',2,'verbose',0);
t_fm2 = toc;
mse_fm2 = (norm(Wgen - Result2.x,'fro')/norm(Wgen,'fro'))^2;
%=== BSBL-FM ==============================================================

fprintf('BSBL-FM(learn correlation) : time: %4.3f, MSE: %g, Iter=%d\n',mean(t_fm2),mean(mse_fm2),Result2.count);

%=== draw(1)
figure(1)
clf;
subplot(121)
plot(Wgen,'b-','linewidth',2); hold on; grid on; axis tight
hx1 = xlabel('(a) Original Signal'); hy1 = ylabel('Amplitude');
ax1 = gca;
subplot(122)
plot(Result2.x,'b-','linewidth',2); hold on; grid on; axis tight
hx2 = xlabel('(b) Recover by BSBL-FM');
ax2 = gca;

%--- config ---
set(ax1, 'LooseInset', get(ax1, 'TightInset'));
set(ax2, 'LooseInset', get(ax2, 'TightInset'));
set([ax1 ax2 hl],'FontName','Times','FontSize',13);
set([hx1 hy1 hx2],'FontName','Times','FontSize',15,'FontWeight','bold');

