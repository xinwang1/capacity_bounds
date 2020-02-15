% This script computes the upper bound in Eq.(8) of arXiv:1912.00931.
% Required packages: 
% CVX http://cvxr.com/cvx/download/
% QETLAB http://www.qetlab.com/Main_Page
% CVX_quad https://github.com/hfawzi/cvxquad

%% input the depolarizing channel
p = 0.05/(3/4); %noise parameter for the channel D_p(\rho):=(1-p)\rho + p tr(\rho)I_d/2
JD = (1 - p)*MaxEntangled(2,0,1)*MaxEntangled(2,0,1)' + p*eye(4)/4; %choi matrix of the channel

s_min = 1;  % assume that Q1 is 1
e_c = 0.02; % the accuracy of the brute-force search
for c=0:e_c:1  % brute-force search over the parameter of the extended states
S0 = [c sqrt(c*(1-c));sqrt(c*(1-c)) 1-c]; %flag 1
S1 = [1 0;0 0]; % flag 2
JF = (1 - p)*kron(MaxEntangled(2,0,1)*MaxEntangled(2,0,1)',S0) + p*kron(eye(4)/4,S1); %the choi state of flagged channel with pure flag states
[epsilon,de,JFC] = epsilon_degradable(2,4,JF); %optimize the epsilon degradable parameter
s_info = -quantum_cond_entr(JF,[2 4],1)/log(2) +  4*epsilon*log2(de) ...
            +  2*BosonicEntropy(epsilon); ;  %bound via the epsilon degradable state
if  s_info < s_min 
    cmin = c;
    s_min = s_info;  %choose the smaller one
    JFC_min = JFC;
else
end  
end

Q_bound = s_min


%% function of binary entropy
function entropy = BinaryEntropy(lambda)
entropy = - lambda*log2(lambda) - (1-lambda)*log2((1-lambda));
end
%% function of bosonic entropy
function entropy = BosonicEntropy(lambda)
entropy = (1+lambda)*BinaryEntropy(lambda/(1+lambda));
end

  %% compute the epsilon degradability
function [eps,de,J] = epsilon_degradability(da,db,JN)
JC = ComplementaryMap(JN,[da db]); %the complementary channel
rhoAE = JC/da;
rhoAB = JN/da;
s = size(JC);
de = s(1)/da; %dimension of environment
cvx_begin sdp 
    variable J(db*de,db*de) hermitian;
    sigmaAE = PartialTrace(kron(eye(da),J)*kron(PartialTranspose(rhoAB,2,[da,db]),eye(de)),2,[da,db,de]);
    eps = TraceNorm(sigmaAE - rhoAE)/2;
    minimize eps
    subject to
            J >= 0;
            PartialTrace(J,2,[db,de]) == eye(db);
   cvx_end
end
