function Result = BSBL_BO(Phi, y, blkStartLoc, LearnLambda, varargin)

% BSBL-BO: Recover block sparse signal (1D) exploiting intra-block correlation, given the block partition.
%
%          The algorithm solves the inverse problem for the block sparse
%          model with known block partition:
%                        y = Phi * x + v
%
%
% ============================== INPUTS ============================== 
%   Phi         : N X M known matrix
%
%   y           : N X 1 measurement vector 
%
%   blkStartLoc : Start location of each block
%   
%   LearnLambda : (1) If LearnLambda = 1, use the lambda learning rule for very LOW SNR cases (SNR<10dB)
%                     (using lambda=std(y)*1e-2 or user-input value as initial value)
%                 (2) If LearnLambda = 2, use the lambda learning rule for medium noisy cases (SNR>10dB) 
%                     (using lambda=std(y)*1e-2 or user-input value as initial value)
%                 (3) If LearnLambda = 0, do not use the lambda learning rule 
%                     ((using lambda=1e-14 or user-input value as initial value)
%                 
%
% [varargin values -- in most cases you can use the default values]
%
%   'LEARNTYPE'  : LEARNTYPE = 0: Ignore intra-block correlation
%                  LEARNTYPE = 1: Exploit intra-block correlation 
%                 [ Default: LEARNTYPE = 1 ]
%
%  'PRUNE_GAMMA'  : threshold to prune out small gamma_i 
%                   (generally, 10^{-3} or 10^{-2})
%
%  'LAMBDA'       : user-input value for lambda
%                  [ Default: LAMBDA=1e-14 when LearnLambda=0; LAMBDA=std(y)*1e-2 in noisy cases]
%
%  'MAX_ITERS'    : Maximum number of iterations.
%                 [ Default value: MAX_ITERS = 600 ]
%
%  'EPSILON'      : Solution accurancy tolerance parameter 
%                 [ Default value: EPSILON = 1e-8   ]
%
%  'PRINT'        : Display flag. If = 1: show output; If = 0: supress output
%                 [ Default value: PRINT = 0        ]
%
% ==============================  OUTPUTS ============================== 
%   Result : 
%      Result.x          : the estimated block sparse signal
%      Result.gamma_used : indexes of nonzero groups in the sparse signal
%      Result.gamma_est  : the gamma values of all the groups of the signal
%      Result.B          : the final value of the B
%      Result.count      : iteration times
%      Result.lambda     : the final value of lambda
%
%
% ========================= Command examples  =============================
%   < Often-used command >
%    For most noisy environment (SNR > 10dB):
%          
%           Result =  BSBL_BO(Phi, y, blkStartLoc, 2);  
%
%    For very low SNR cases (SNR < 10 dB):
%           
%           Result =  BSBL_BO(Phi, y, blkStartLoc, 1);   
%
%    For noiseless cases:
%          
%           Result =  BSBL_BO(Phi, y, blkStartLoc, 0);  
%
%    To recover non-Sparse structured signals (noiseless):
%           Result = BSBL_BO(Phi,y,groupStartLoc,0,'prune_gamma',-1);
%                    ('prune_gamma' can be set any non positive constant)
%
%    < Full-Command Example >
%           Result =  BSBL_BO(Phi, y, blkStartLoc, learnlambda, ...
%                                                 'LEARNTYPE', 1,...
%                                                 'PRUNE_GAMMA',1e-2,...
%                                                 'LAMBDA',1e-3,...
%                                                 'MAX_ITERS', 800,...
%                                                 'EPSILON', 1e-8,...
%                                                 'PRINT',0);
%
% ================================= See Also =============================
%   EBSBL_BO,   BSBL_EM,  BSBL_L1,  EBSBL_L1,  TMSBL,    TSBL      
%
% ================================ Reference =============================
%   [1] Zhilin Zhang, Bhaskar D. Rao, Extension of SBL Algorithms for the 
%       Recovery of Block Sparse Signals with Intra-Block Correlation, 
%       available at: http://arxiv.org/abs/1201.0862
%  
%   [2] Zhilin Zhang, Tzyy-Ping Jung, Scott Makeig, Bhaskar D. Rao, 
%       Low Energy Wireless Body-Area Networks for Fetal ECG Telemonitoring 
%       via the Framework of Block Sparse Bayesian Learning, 
%       available at: http://arxiv.org/pdf/1205.1287v1.pdf
%
%   [3] webpage: http://dsp.ucsd.edu/~zhilin/BSBL.html, or
%                https://sites.google.com/site/researchbyzhang/bsbl
%
% ============= Author =============
%   Zhilin Zhang (z4zhang@ucsd.edu, zhangzlacademy@gmail.com)
%
% ============= Version =============
%   1.4 (07/23/2012) debug
%   1.3 (05/30/2012) make faster
%   1.2 (05/28/2012)
%   1.1 (01/22/2012)
%   1.0 (08/27/2011)
%


% scaling...
scl = std(y);
if (scl < 0.4) | (scl > 1)
    y = y/scl*0.4;
end

% Default Parameter Values for Any Cases
EPSILON       = 1e-8;       % solution accurancy tolerance
MAX_ITERS     = 600;        % maximum iterations
PRINT         = 0;          % don't show progress information
LEARNTYPE     = 1;          % adaptively estimate the covariance matrix B

if LearnLambda == 0  
    lambda = 1e-12;   
    PRUNE_GAMMA = 1e-3;
elseif LearnLambda == 2
    lambda = scl * 1e-2;    
    PRUNE_GAMMA = 1e-2;
elseif LearnLambda == 1
    lambda = scl * 1e-2;    
    PRUNE_GAMMA = 1e-2;
else
    error(['Unrecognized Value for Input Argument ''LearnLambda''']);
end


if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'learntype'
                LEARNTYPE = varargin{i+1};
                if LEARNTYPE ~= 1 & LEARNTYPE ~= 0
                    error(['Unrecognized Value for Input Argument ''LEARNTYPE''']);
                end
            case 'prune_gamma'
                PRUNE_GAMMA = varargin{i+1}; 
            case 'lambda'
                lambda = varargin{i+1};    
            case 'epsilon'   
                EPSILON = varargin{i+1}; 
            case 'print'    
                PRINT = varargin{i+1}; 
            case 'max_iters'
                MAX_ITERS = varargin{i+1};  
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end


if PRINT
    fprintf('\n====================================================\n');
    fprintf('           Running BSBL-BO ....... \n');
    fprintf('           Information about parameters...\n');
    fprintf('====================================================\n');
    fprintf('PRUNE_GAMMA  : %e\n',PRUNE_GAMMA);
    fprintf('lambda       : %e\n',lambda);
    fprintf('LearnLambda  : %d\n',LearnLambda);    
    fprintf('LearnType    : %d\n',LEARNTYPE);
    fprintf('EPSILON      : %e\n',EPSILON);
    fprintf('MAX_ITERS    : %d\n\n',MAX_ITERS);
end


%% Initialization
[N,M] = size(Phi);
Phi0 = Phi;
blkStartLoc0 = blkStartLoc;
p = length(blkStartLoc);   % block number
for k = 1 : p-1
    blkLenList(k) = blkStartLoc(k+1)-blkStartLoc(k);
end
blkLenList(p) = M - blkStartLoc(end)+1;
maxLen = max(blkLenList);
if sum(blkLenList == maxLen) == p, 
    equalSize = 1;
else
    equalSize = 0;
end

for k = 1 : p
    Sigma0{k} = eye(blkLenList(k));
end

gamma = ones(p,1);
keep_list = [1:p]';
usedNum = length(keep_list);
mu_x = zeros(M,1);
count = 0;


%% Iteration
while (1)
    count = count + 1;

    %=========== Prune weighys as yheir hyperparameyers go yo zero ==============
    if (min(gamma) < PRUNE_GAMMA)
        index = find(gamma > PRUNE_GAMMA);
        usedNum = length(index);
        keep_list = keep_list(index); 
        if isempty(keep_list), 
            fprintf('\n====================================================================================\n');
            fprintf('x becomes zero vector. The solution may be incorrect. \n');
            fprintf('Current ''prune_gamma'' = %g, and Current ''EPSILON'' = %g.\n',PRUNE_GAMMA,EPSILON);
            fprintf('Try smaller values of ''prune_gamma'' and ''EPSILON'' or normalize ''y'' to unit norm.\n');
            fprintf('====================================================================================\n\n');
            break; 
        end;
        blkStartLoc = blkStartLoc(index);
        blkLenList = blkLenList(index);
        
        % prune gamma and associated components in Sigma0 
        gamma = gamma(index);
        temp = Sigma0;
        Sigma0 = [];
        for k = 1 : usedNum
            Sigma0{k} = temp{index(k)};
        end
        
        % construct new Phi
        temp = [];
        for k = 1 : usedNum
            temp = [temp, Phi0(:,blkStartLoc(k):blkStartLoc(k)+blkLenList(k)-1)];
        end
        Phi = temp;
        %clear temp;
    end

    %=================== Compute new weights =================
    mu_old = mu_x;
    
    PhiBPhi = zeros(N);
    currentLoc = 0;
    for i = 1 : usedNum
        
        currentLen = size(Sigma0{i},1);
        currentLoc = currentLoc + 1;
        currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
        
        PhiBPhi = PhiBPhi + Phi(:, currentSeg)*Sigma0{i}*Phi(:, currentSeg)';
        currentLoc = currentSeg(end);
    end
     
    H = Phi' /(PhiBPhi + lambda * eye(N));
    Hy = H * y;      
    HPhi = H * Phi;
    
    mu_x = zeros(size(Phi,2),1);
    Sigma_x = [];
    Cov_x = [];
     
    B = []; invB = []; B0 = zeros(maxLen); r0 = zeros(1); r1 = zeros(1);
    currentLoc = 0;
    for i = 1 : usedNum
        
        currentLen = size(Sigma0{i},1);
        currentLoc = currentLoc + 1;
        seg = currentLoc : 1 : currentLoc + currentLen - 1;
        
        mu_x(seg) = Sigma0{i} * Hy(seg);       % solution
        Sigma_x{i} = Sigma0{i} - Sigma0{i} * HPhi(seg,seg) * Sigma0{i};
        Cov_x{i} = Sigma_x{i} + mu_x(seg) * mu_x(seg)';
        currentLoc = seg(end);
        
        %=========== Learn correlation structure in blocks ===========
        % do not consider correlation structure in each block
        if LEARNTYPE == 0
            B{i} = eye(currentLen);
            invB{i} = eye(currentLen);
 
        % constrain all the blocks have the same correlation structure
        elseif LEARNTYPE == 1
            if equalSize == 0
                if currentLen > 1
                    temp = Cov_x{i}/gamma(i);
                    r0 = r0 + mean(diag(temp));
                    r1 = r1 + mean(diag(temp,1));
                end
            elseif equalSize == 1
                B0 = B0 + Cov_x{i}/gamma(i);
            end

        end % end of learnType

    end
    
    %=========== Learn correlation structure in blocks with Constraint 1 ===========
    % If blocks have the same size
    if (equalSize == 1) & (LEARNTYPE == 1)

        % Constrain all the blocks have the same correlation structure
        % (an effective strategy to avoid overfitting)
        b = (mean(diag(B0,1))/mean(diag(B0)));
        if abs(b) >= 0.99, b = 0.99*sign(b); end;
        bs = [];
        for j = 1 : maxLen, bs(j) = (b)^(j-1); end;
        B0 = toeplitz(bs);
 
        for i = 1 : usedNum
             
            B{i} = B0;
            invB{i} = inv(B{i});
        end
    
    % if blocks have different sizes
    elseif (equalSize == 0) & (LEARNTYPE == 1)
        r = r1/r0; if abs(r) >= 0.99, r = 0.99*sign(r); end;

        for i = 1 : usedNum
            currentLen = size(Cov_x{i},1);

            bs = [];
            for j = 1 : currentLen, bs(j) = r^(j-1); end;
            B{i} = toeplitz(bs);
            invB{i} = inv(B{i});

        end
         
    end

    
    % estimate gamma(i) and lambda 
    if LearnLambda == 1          
        gamma_old = gamma;
        lambdaComp = 0; currentLoc = 0;
        for i =  1 : usedNum

            currentLen = size(Sigma_x{i},1);
            currentLoc = currentLoc + 1;
            currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
            
            gamma(i) = gamma_old(i)*norm( sqrtm(B{i})*Hy(currentSeg) )/sqrt(trace(HPhi(currentSeg,currentSeg)*B{i}));
            
            lambdaComp = lambdaComp + trace(Phi(:,currentSeg)*Sigma_x{i}*Phi(:,currentSeg)');
            
            Sigma0{i} = B{i} * gamma(i);
            
            currentLoc = currentSeg(end);
        end
        lambda = norm(y - Phi * mu_x,2)^2/N + lambdaComp/N; 
        
          
    elseif LearnLambda == 2    
        gamma_old = gamma;
        lambdaComp = 0;  currentLoc = 0;
        for i =  1 : usedNum  
            
            currentLen = size(Sigma_x{i},1);
            currentLoc = currentLoc + 1;
            currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
 
            gamma(i) = gamma_old(i)*norm( sqrtm(B{i})*Hy(currentSeg) )/sqrt(trace(HPhi(currentSeg,currentSeg)*B{i}));
            
            lambdaComp = lambdaComp + trace(Sigma_x{i}*invB{i})/gamma_old(i);
            
            Sigma0{i} = B{i} * gamma(i);
            
            currentLoc = currentSeg(end);
        end
        lambda = norm(y - Phi * mu_x,2)^2/N + lambda * (length(mu_x) - lambdaComp)/N; 
        
    else   % only estimate gamma(i)
        gamma_old = gamma;
        currentLoc = 0;
        for i =  1 : usedNum
            % gamma(i) = trace(invB{i} * Cov_x{i})/size(Cov_x{i},1);
            currentLen = size(Sigma0{i},1);
            currentLoc = currentLoc + 1;
            seg = currentLoc : 1 : currentLoc + currentLen - 1;
            
            gamma(i) = gamma_old(i)*norm( sqrtm(B{i})*Hy(seg) )/sqrt(trace(HPhi(seg,seg)*B{i}));
            
            Sigma0{i} = B{i} * gamma(i);
            
            currentLoc = seg(end);
        end
    end

    
    % ================= Check stopping conditions, eyc. ==============
    if (size(mu_x) == size(mu_old))
        dmu = max(max(abs(mu_old - mu_x)));
        if (dmu < EPSILON)  break;  end;
    end;
    if (PRINT) 
        disp([' iters: ',num2str(count),...
            ' num coeffs: ',num2str(usedNum), ...
            ' min gamma: ', num2str(min(gamma)),...
            ' gamma change: ',num2str(max(abs(gamma - gamma_old))),...
            ' mu change: ', num2str(dmu)]); 
    end;
    if (count >= MAX_ITERS), if PRINT, fprintf('Reach max iterations. Stop\n\n'); end; break;  end;

end;



if isempty(keep_list)
    Result.x = zeros(M,1);
    Result.gamma_used = [];
    Result.gamma_est = zeros(p,1);
    Result.B = B;
    Result.count = count;
    Result.lambdatrace = lambda;

else
    %% Expand hyperparameyers
    gamma_used = sort(keep_list);
    gamma_est = zeros(p,1);
    gamma_est(keep_list,1) = gamma;


    %% reconstruct the original signal
    x = zeros(M,1);
    currentLoc = 0;
    for i = 1 : usedNum

        currentLen = size(Sigma0{i},1);
        currentLoc = currentLoc + 1;
        seg = currentLoc : 1 : currentLoc + currentLen - 1;

        realLocs = blkStartLoc0(keep_list(i)) : blkStartLoc0(keep_list(i))+currentLen-1;

        x( realLocs ) = mu_x( seg );
        currentLoc = seg(end);
    end

    if (scl < 0.4) | (scl > 1)
        Result.x = x * scl/0.4;
    else
        Result.x = x;
    end
    Result.gamma_used = gamma_used;
    Result.gamma_est = gamma_est;
    Result.B = B;
    Result.count = count;
    Result.lambda = lambda;
end

return;



