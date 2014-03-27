% MMV for BSBL-FM
clear;  close all;

% get the raw dataset (which includes strong noise)
load signal_01.mat;

% downsampling to 250 Hz
s = s(:,1:4:2048); 
C = 8;
 
% the size of the sensing matrix Phi
load Phi;
[M,N] = size(Phi);
W = dctmtx(N); W = W';
A = Phi*W;
% A=zeros(M,N);
% for k=1:M
%     A(k,:)=dct(Phi(k,:));
% end

% block size of the user-defined block partition of BSBL-BO
blkLen = 32;

% the user-defined block partition
groupStartLoc = 1:blkLen:N;

% variable for recovered dataset
sig = s';

%====================================================
% Compress all the ECG recordings, and recover them
%====================================================
Y = Phi*sig;

% method 1. recover using MBSBL-FM
tic;
Result = BSBL_FM(A,Y,groupStartLoc,99,'epsilon',1e-4,'learnType',0,'verbose',0);
runtime = toc;

% method 2. recover using SMV each channel
runtime1 = zeros(C,1); mse1 = zeros(C,1);
for ii = 1 : C
    y = Phi*sig(:,ii);
    tic;
    Result1 = BSBL_FM(A,y,groupStartLoc,99,'epsilon',1e-4,'learnType',0,'verbose',0);
    runtime1(ii) = toc;
    
    mse1(ii) = (norm(sig(:,ii) - W*Result1.x,'fro')/norm(sig(:,ii),'fro'))^2;
end

fprintf('Tol Runtime %8.8f <----> Single Tol Runtime %8.8f \n',runtime,sum(runtime1));
for jj = 1 : 8 % 8 channels
    mse = (norm(sig(:,jj) - W*Result.x(:,jj),'fro')/norm(sig(:,jj),'fro'))^2;
    fprintf('MSE(C=%3d): %8.8f <----> Single MSE: %8.8f\n',jj,mse,mse1(jj)); 
end
