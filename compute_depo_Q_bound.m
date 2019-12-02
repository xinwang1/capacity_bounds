% This script computes the upper bound Q_pf in the paper.
% Required package: 
% CVX http://cvxr.com/cvx/download/
% QETLAB http://www.qetlab.com/Main_Page
% CVX_quad https://github.com/hfawzi/cvxquad

%% input the depolarizing channel
p = 0.05/(3/4); %noise parameter for the channel D_p(\rho):=(1-p)\rho + p tr(\rho)I_d/2
JD = (1 - p)*MaxEntangled(2,0,1)*MaxEntangled(2,0,1)' + p*eye(4)/4; %choi matrix of the channel

IC_min = 1;  % assume that Q1 is 1
e_c = 0.02; %the accuracy of the brute-force search
for c=0:e_c:1
S0 = [c sqrt(c*(1-c));sqrt(c*(1-c)) 1-c]; %flag 1
S1 = [1 0;0 0]; % flag 2
JF = (1 - p)*kron(MaxEntangled(2,0,1)*MaxEntangled(2,0,1)',S0) + p*kron(eye(4)/4,S1); %the choi matrix of flagged channel with pure flag states
[epsilon,de,JFC] = epsilon_degradable(2,4,JF); %optimize the epsilon degradable channel

coh_info = -quantum_cond_entr(JF,[2 4],1)/log(2) + epsilon*log2(de-1) + 2*epsilon*log2(de) ...
            + BinaryEntropy(epsilon) + BosonicEntropy(epsilon);  %bound via the epsilon degradable channel
        
if  coh_info < IC_min 
    cmin = c;
    IC_min = coh_info;  %choose the smaller one
    JFC_min = JFC;
else
end  
end

Q_pf = IC_min
%% function of binary entropy
function entropy = BinaryEntropy(lambda)
entropy = - lambda*log2(lambda) - (1-lambda)*log2((1-lambda));
end
%% function of bosonic entropy
function entropy = BosonicEntropy(lambda)
entropy = (1+lambda)*BinaryEntropy(lambda/(1+lambda));
end
%% optimize the epsilon degradable channel
function [eps,de,J] = epsilon_degradable(da,db,JN)
JC = ComplementaryMap(JN,[da db]); %the complementary channel
s = size(JC);
de = s(1)/da; %dimension of environment
cvx_begin sdp 
    variable Z(da*de,da*de) hermitian;
    variable J(db*de,db*de) hermitian;
    Jcompose = PartialTrace(kron(eye(da),J)*kron(PartialTranspose(JN,2,[da,db]),eye(de)),2,[da,db,de]);
    eps = SchattenNorm(PartialTrace(Z,2,[da de]),inf);
    minimize eps
    subject to
            Z >= JC - Jcompose;
            Z >= 0;
            J >= 0;
            PartialTrace(J,2,[db,de]) == eye(db);
   cvx_end
end

