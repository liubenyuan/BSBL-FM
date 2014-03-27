function Result = BSBL_FM(PHI,y,blkStartLoc,LearnLambda,varargin)
%------------------------------------------------------------------
% The block BCS algorithm for our following paper:
%   "Fast Marginalized Block SBL Algorithm" (Preprint, 2012)
%
% for Zhang Zhilin's 
%   "Extension of SBL Algorithms for the Recovery of Block 
%    Sparse Signals with Intra-Block Correlation" (Preprint, Zhang2012)
%
% Coded by: Liu Benyuan
% Change Log:
%     v1.5[20121122]: optimized for speed
%     v1.6[20121122]: add complex support, only works for learnType=0;
%     v1.7[20121126]: add comments
%
%------------------------------------------------------------------
% Input for BSBL-FM:
%   PHI: projection matrix
%   y:   CS measurements
%   blkStartLoc : Start location of each block
%   LearnLambda : (1) If LearnLambda = 1, 
%                     use the lambda learning rule for MEDIUM SNR cases (SNR<=30dB)
%                     (using lambda=std(y)*1e-2 or user-input value as initial value)
%                 (2) If LearnLambda = 2, 
%                     use the lambda learning rule for HIGH SNR cases (SNR>30dB) 
%                     (using lambda=std(y)*1e-3 or user-input value as initial value)
%                 (3) If LearnLambda = 0, do not use the lambda learning rule 
%                     ((using lambda=1e-7 or user-input value as initial value)
%
% [varargin values -- in most cases you can use the default values]
%   'LEARNTYPE'  : LEARNTYPE = 0: Ignore intra-block correlation
%                  LEARNTYPE = 1: Exploit intra-block correlation 
%                 [ Default: LEARNTYPE = 1 ]
%   'VERBOSE'    : debuging information.
%   'EPSILON'    : convergence criterion
%
% ==============================  OUTPUTS ============================== 
%   Result : 
%      Result.x          : the estimated block sparse signal
%      Result.gamma_used : indexes of nonzero groups in the sparse signal
%      Result.gamma_est  : the gamma values of all the groups of the signal
%      Result.B          : the final mean value of each correlation block
%      Result.count      : iteration times
%      Result.lambda     : the final value of lambda

% default values for BSBL-FM
eta = 1e-4;      % default convergence test
verbose = 0;     % print some debug information
learnType = 0;   % default not to exploit intra block correlation
max_it = 1000;   % maximum iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. intialize, scale
scl = max(std(y)); % max scale
if (scl < 0.4) || (scl > 1)
    y = y/scl*0.4;
end
[~,M] = size(PHI);
[~,T] = size(y);

% select sigma2
stdy2 = mean(std(y))^2;
sigma2 = 1e-3*stdy2;       % default value if otherwise specified [99]
if LearnLambda == 0 
	sigma2 = 1e-6;         % noiseless                            [0 ]
elseif LearnLambda == 2
	sigma2 = 1e-2*stdy2;   % high SNR (SNR>=20)                   [2 ]
elseif LearnLambda == 1
	sigma2 = 1e-1*stdy2;   % medium SNR (SNR<20)                  [1 ]
end

if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'learntype'
				learnType = varargin{i+1};
			case 'epsilon'
				eta = varargin{i+1};
            case 'sigma2_scale'
                sigma2 = varargin{i+1}*stdy2;
            case 'max_iters'
                max_it = varargin{i+1};
            case 'verbose'    
                verbose = varargin{i+1}; 
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end

beta = 1/sigma2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. formalize the blocks and quantities used in the code
%    p           : the number of blocks
%    blkStartLoc : the start index of blk
%    blkLenList  : the length of each block
p = length(blkStartLoc);
blkLenList = ones(p,1);
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
% when the blkLen=1 we avoid the exploiting feature.
if maxLen == 1,
	learnType = 0;
end

% pre-allocating space
S          = cell(p,1); s = cell(p,1);
Q          = cell(p,1); q = cell(p,1);
currentSeg = cell(p,1);
localSeg   = cell(p,1);
Phi        = cell(p,1);
% 2. prepare the quantities used in the code.
for k = 1 : p
	currentLoc    = blkStartLoc(k);
	currentLen    = blkLenList(k);
	currentSeg{k} = currentLoc:1:currentLoc + currentLen - 1;

	Phi{k} = PHI(:,currentSeg{k});
	S{k}   = beta.*Phi{k}'*Phi{k};
	Q{k}   = beta.*Phi{k}'*y;
end

% 3. start from *NULL*, decide which one to add ->
A     = cell(p,1); 
Am    = cell(p,1); % old A
theta = zeros(p,1);
for k = 1 : p
	A{k}     = (S{k})\(Q{k}*Q{k}' - S{k})/(S{k});
	theta(k) = 1/blkLenList(k) * real(trace(A{k}));
	A{k}     = eye(blkLenList(k)).*theta(k);
end
% select the basis that minimize the change of *likelihood*
ml  = inf*ones(1,p);
ig0 = find(theta>0);
len = length(ig0);
for kk = 1:len
	k = ig0(kk);
	ml(k) = log(abs(det(eye(blkLenList(k)) + A{k}*S{k}))) ...
	      - trace(real(Q{k}'/(eye(blkLenList(k)) + A{k}*S{k})*A{k}*Q{k}));
end
[~,index] = min(ml);
gamma = theta(index);
Am{index} = A{index}; % Am -> record the past value of A
if verbose, fprintf(1,'ADD,\t idx=%3d, GAMMA_OP=%f\n',index,gamma); end

% 3. update quantities (Sig,Mu,S,Q,Phiu) 
Sigma_ii = (eye(blkLenList(index))/Am{index} + S{index})\eye(blkLenList(index));
Sig      = Sigma_ii;
Mu       = Sigma_ii*Q{index};
% The relevent block basis
Phiu = Phi{index};
for k = 1 : p
	Phi_k = Phi{k};
	S{k}  = S{k} - beta^2.*Phi_k'*(Phiu*Sigma_ii*Phiu')*Phi_k;
	Q{k}  = Q{k} - beta  .*Phi_k'*Phiu*Mu;
end

% system parameter
ML=zeros(max_it,1);

for count = 1:max_it
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	localLoc = 1;
	for i = 1 : length(index);
		k = index(i);
		localLen = blkLenList(k);
		localSeg{i} = localLoc:1:localLoc + localLen - 1;
		localLoc = localLoc + localLen;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% pre-process steps if we want to learn the intra-block-correlation
	%     learnType == 2 : calculate the mean of r_i
	if learnType == 2
		len = length(index); r = zeros(len,1);
		for i = 1 : len
			seg      = localSeg{i};
			Sigma_ii = Sig(seg,seg);
			Mu_i     = Mu(seg);
			[~,r(i)] = learnB(Sigma_ii,Mu_i,gamma(i));
		end
		r_hat = mean(r); % mean or max
		BT    = genB(r_hat,maxLen);
	end

	% calculate s,q
	for k = 1 : p
		which = find(index==k,1);
        if isempty(which)   % the k-th basis is not included
			s{k} = S{k};
			q{k} = Q{k};
        else                % the k-th basis is calculated
			invDenom = (eye(blkLenList(k)) - S{k}*Am{k})\eye(blkLenList(k));
			s{k} = invDenom*S{k};
			q{k} = invDenom*Q{k};
        end
		% learnType ==>> [0,1,2]
 		A{k} = (s{k})\(q{k}*q{k}' - s{k})/(s{k});
		theta(k) = 1/blkLenList(k) * real(trace(A{k}));
		if learnType == 0      % [0] without intra-correlation
			A{k} = eye(blkLenList(k))*theta(k);
		elseif learnType == 1  % [1] with individual intra corr
			rr = mean(diag(A{k},1))/mean(diag(A{k}));
			if abs(rr)>0.95, rr = 0.95*sign(rr); end
			Bc = genB(rr,blkLenList(k));
			A{k} = Bc*theta(k);
		elseif learnType == 2  % [2] with unified intra corr
			if equalSize 
				Bc = BT;
			else
				Bc = genB(r_hat,blkLenList(k));
			end
			A{k} = Bc.*theta(k);
		end
	end

    % choice the next basis that [minimizes] the cost function
    ml =  inf*ones(1,p);
    ig0 = find(theta>0);
    % index for re-estimate
    [ire,~,which] = intersect(ig0,index);
    if ~isempty(ire)
		len = length(which);
		for kk = 1:len
			k = ire(kk);
			ml(k) = log(abs(det(eye(blkLenList(k)) + A{k}*s{k}))) ...
			        -trace(real(q{k}'/(eye(blkLenList(k)) + A{k}*s{k})*A{k}*q{k})) ...
				   -(log(abs(det(eye(blkLenList(k))+ Am{k}*s{k}))) ...
				    -trace(real(q{k}'/(eye(blkLenList(k)) + Am{k}*s{k})*Am{k}*q{k})));
		end
    end
    % index for adding
    iad = setdiff(ig0,ire);
    if ~isempty(iad)
		len = length(iad);
		for kk = 1:len
			k = iad(kk);
			ml(k) = log(abs(det(eye(blkLenList(k)) + A{k}*s{k}))) ...
				   -trace(real(q{k}'/(eye(blkLenList(k)) + A{k}*s{k})*A{k}*q{k}));
		end
    end
    % index for deleting
    is0 = setdiff((1:p),ig0);
    [ide,~,which] = intersect(is0,index);
    if ~isempty(ide)
		len = length(which);
		for kk = 1:len
			k = ide(kk);
			ml(k) = -(log(abs(det(eye(blkLenList(k)) + Am{k}*s{k}))) ...
			         -trace(real(q{k}'/(eye(blkLenList(k)) + Am{k}*s{k})*Am{k}*q{k})));
		end
    end

	% as we are minimizing the cost function :
    [ML(count),idx] = min(ml);
    
    % check if terminates?
	if ML(count)>=0, break; end
    if count > 2 && abs(ML(count)-ML(count-1)) < abs(ML(count)-ML(1))*eta, break; end

    % update block gammas
    which = find(index==idx);
	% processing the quantities update
	if ~isempty(which)  % the select basis is already in the *LIST*
		seg    = localSeg{which};
		Sig_j  = Sig(:,seg);
		Sig_jj = Sig(seg,seg);
		if theta(idx)>0
			%%%% re-estimate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if verbose,fprintf(1,'REE,\t idx=%3d, GAMMA_OP=%f\n',idx,theta(idx));end
			gamma_new = theta(idx);
			ki  = Sig_j/(Sig_jj + Am{idx}/(Am{idx} - A{idx})*A{idx})*Sig_j';
			Sig = Sig - ki;
			Mu  = Mu - beta.*ki*Phiu'*y;
			PKP = Phiu*ki*Phiu';
            for k = 1 : p
				Phi_m = Phi{k};
				PPKP  = Phi_m'*PKP;
				S{k}  = S{k} + beta^2.*PPKP*Phi_m;
				Q{k}  = Q{k} + beta^2.*PPKP*y;
            end
			%
			gamma(which) = gamma_new; % 1
			Am{idx}      = A{idx};    % 2
		else 
			%%%% delete %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if verbose,fprintf(1,'DEL,\t idx=%3d, GAMMA_OP=%f\n',idx,gamma(which));end
			if length(index)==1, break; end % we are deleting the only one
			ki  = Sig_j/Sig_jj*Sig_j';
			Sig = Sig - ki;
			Mu  = Mu - beta.*ki*Phiu'*y;
			PKP = Phiu*ki*Phiu';
			for k = 1 : p
				Phi_m = Phi{k};
				PPKP  = Phi_m'*PKP;
				S{k}  = S{k} + beta^2.*PPKP*Phi_m;
				Q{k}  = Q{k} + beta^2.*PPKP*y;
			end
			% delete relevant basis and block
            index(which) = [];
			Mu(seg,:)    = [];
			Sig(:,seg)   = [];
			Sig(seg,:)   = [];
			Phiu(:,seg)  = [];
			%
			gamma(which) = []; % 1
			Am{idx}      = []; % 2
		end
	else
		if theta(idx)>0
			%%%% add %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if verbose,fprintf(1,'ADD,\t idx=%3d, GAMMA_OP=%f\n',idx,theta(idx));end
			gamma_new = theta(idx);
			Phi_j     = Phi{idx};
			%
			Sigma_ii = (eye(blkLenList(idx))+A{idx}*S{idx})\A{idx};
			mu_i     = Sigma_ii*Q{idx};
			SPP      = Sig*Phiu'*Phi_j; % common
			Sigma_11 = Sig + beta^2.*SPP*Sigma_ii*SPP';
			Sigma_12 = -beta.*SPP*Sigma_ii;
			Sigma_21 = Sigma_12';
			mu_1     = Mu - beta.*SPP*mu_i;
			e_i      = Phi_j - beta.*Phiu*SPP;
			ESE      = e_i*Sigma_ii*e_i';
			for k = 1 : p
				Phi_m = Phi{k};
				S{k}  = S{k} - beta^2.*Phi_m'*ESE*Phi_m;
				Q{k}  = Q{k} - beta.*Phi_m'*e_i*mu_i;
			end
			% adding relevant basis
			Sig     = [Sigma_11 Sigma_12; ...
			           Sigma_21 Sigma_ii];
			Mu      = [mu_1; ...
			           mu_i];
			Phiu    = [Phiu Phi_j];
			index   = [index;idx];
			gamma   = [gamma;gamma_new];  % 1
			Am{idx} = A{idx};             % 2
		else
			break; % null operation
		end
	end

end
% format the output ===> X the signal
weights = zeros(M,T);
formatSeg = [currentSeg{index}];
weights(formatSeg,:) = Mu;
if (scl < 0.4) || (scl > 1)
    Result.x = weights * scl/0.4;
else
    Result.x = weights;
end
Result.r = 1.0; % lazy ...
Result.gamma_used = index;
Result.gamma_est = gamma;
Result.count = count;
Result.lambda = sigma2;
% END %

%% sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions of estimating the AR(1) coefficient r and 
% reconstruction the covariance matrix with B^{-1} valid
function [B,r] = learnB(Sig,Mu,gamma)
len = length(Mu);
B = (Sig + Mu*Mu')./gamma;
r = (mean(diag(B,1))/mean(diag(B)));
if abs(r) >= 0.95, r = 0.95*sign(r); end;
B = genB(r,len);
% generate B according to r,len 
% NOTE: abs(r) should be less than 1.0
function B = genB(r,len)
jup = 0:len-1;
bs = r.^jup;
B = toeplitz(bs);

% generate temporal Smooth matrix
% NOTE: current does not handle L
function B = temporalSmooth(a,b,~,len)
A1 = b.*eye(len);
A2 = (a*b).*[zeros(1,len-1) 0; eye(len-1), zeros(len-1,1)];
Bc = A1 + A2;
B = Bc*Bc';

