% This script computes the upper bound Q_pf in the paper.
% Required package: 
% CVX http://cvxr.com/cvx/download/
% QETLAB http://www.qetlab.com/Main_Page
% CVX_quad https://github.com/hfawzi/cvxquad
% The codes for computing the coherent information and the Rains information of the GAD channel are from  Khatri, S., Sharma, K., & Wilde, M. M. (2019). http://arxiv.org/abs/1903.07747
%% input the depolarizing channel
y = 0.3; N = 0.5; %noise parameter for the GAD channel
A{1}=sqrt(1-N)*[1 0;0 sqrt(1-y)];
A{2}=sqrt(1-N)*[0 sqrt(y);0 0];
A{3}=sqrt(N)*[sqrt(1-y) 0;0 1];
A{4}=sqrt(N)*[0 0;sqrt(y) 0]; %kraus operators of the GAD channel
MES = MaxEntangled(2,0,1)*MaxEntangled(2,0,1)';%normalize max entangled state 
J1 = kron(eye(2),A{1})*MES*kron(eye(2),A{1}');
J2 = kron(eye(2),A{2})*MES*kron(eye(2),A{2}');
J3 = kron(eye(2),A{3})*MES*kron(eye(2),A{3}');
J4 = kron(eye(2),A{4})*MES*kron(eye(2),A{4}');
JA = J1+J2+J3+J4;

S0 = [0 0;0 1]; %flag state 1
S1 = [1 0;0 0]; %flag state 2
JF = kron(J1+J2,S0) + kron(J3+J4,S1); %flagged channel

Q1 = Q1S(JA*2,2,2) 
Q_flagged_channel = Q1S(JF*2,2,4) 
Rains_info = Rains_GADC(y,N)
%% function of lower bound
function value  = Q1S(JN,da,db)
IC_min = 1;
for c=0:0.001:1
    rho = [c 0;0 1-c];
    coh_info = quantum_cond_entr(kron(rho^0.5,eye(db))*JN*kron(rho^0.5,eye(db)),[da db],1)/log(2);
if  coh_info < IC_min 
    IC_min = coh_info;
else
end  
end
value = -IC_min
end

 
 %% function of rains bound
function t  = Rains_GADC(g,N)
% This function is from  http://arxiv.org/abs/1903.07747
t = 0;
for p=0:0.025:1
    Choi_p = [(1-p)*(1-N*g) 0 0 sqrt((1-p)*p*(1-g));
                0 g*N*(1-p) 0 0 ;
                0 0 (1-N)*p*g 0;
                sqrt((1-p)*p*(1-g)) 0 0 p - (1-N)*p*g];
    rains = REE_qubit_X(Choi_p);
if  rains > t 
    t = rains;
else
end  
end
end

function [out] = REE_qubit_X(rho_AB,varargin)
% This function is from  http://arxiv.org/abs/1903.07747

% Last Modified: 14 March 2019

% Calculates the relative entropy of entanglement of a bipartite qubit
% state.

cvx_begin sdp

%cvx_solver sedumi

cvx_expert true

%if ~isempty(varargin)
	cvx_quiet true
%end

variable sigma_AB(4,4) hermitian
variables a b c d  
variable xi complex

RE = quantum_rel_entr(rho_AB,sigma_AB)/log(2);

minimize RE

sigma_AB == 1/2 * [a 0 0 xi; 0 b 0 0; 0 0 c 0; conj(xi) 0 0 d];

sigma_AB >= 0;

trace(sigma_AB) == 1;

PartialTranspose(sigma_AB,2,[2,2]) >= 0;

cvx_end

out=full(RE);

opt_state=full(sigma_AB);

end