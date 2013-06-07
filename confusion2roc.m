function [TPR,FPR,K_1_0,K_0_0, CohensKappa,PPV,NPV, EFF,chi_square,P, Q] = confusion2roc(confusion, sample_size)
% function [TPR,FPR,K_1_0,K_0_0, CohensKappa,PPV,NPV] = confusion2roc(confusion)
%takes a Nx4 confusion matrix of the form [TP,FN,FP,TN] whee these are in
%percentages from 0..1 - i
%sample_size = the number of samples used in creating the confusion matrix
% PPV = Positive predicitive Value
% NPV = negative predictive Value
% EFF = Efficiency
% P = Prevalence
% Q = quality

%the values of confusion need to first be normalized by the sum if any are above 1, where it
%will be assumed that the values represented counts as opposed to
%percentages...

TP = confusion(:,1);
FN = confusion(:,2);
FP = confusion(:,3);
TN = confusion(:,4);

P = TP + FN;
Q = TP + FP;

SE = TP./P;
SP = TN./(1-P);
TPR = SE;
FPR = 1 - SP;

K_1_0 = (SE-Q)./(1-Q);
K_0_0 = (SP-(1-Q))./Q;

nans = isnan(K_1_0(:))|isnan(K_0_0(:));

K_0_0(nans) = 0;

EFF = TP+FN;

if(nargin>1)
    chi_square = sample_size*K_1_0.*K_0_0;
else
    chi_square = {};
end

%page 116 in Kramer's book  Evaluating Medical Tests
% [PQ'*K(1,0)+P'Q*K(0,0)]/(PQ'+P'Q)

CohensKappa = {};
PPV = {};
NPV = {};

if(nargout>=5)
    CohensKappa = (P.*(1-Q).*K_1_0+(1-P).*Q.*K_0_0)./(P.*(1-Q)+(1-P).*Q);
    if(nargout>=6)
        PPV = TP./Q;
        if(nargout>=7)
            NPV = TN./(1-Q);
        end
    end
end
    

