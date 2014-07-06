function [Pfa,Pd,AUC,AUClg,score] = ROC_Curve(tsTgt,tsFa)
% CALL:  [Pfa,Pd,AUC,AUClg,score] = ROC_Curve(tsTgt,tsFa);
%
% PURPOSE: Given two vectors that contain a detector test statistic for
% targets (tsTgt) and false alarms (tsFa), this routine calculates a
% receiver operating characteristic (ROC) curve. To plot a ROC curve, use
% the following commands:
%
%  >> [pfa,pd] = ROC_Curve(tsTgt,tsFa);
%  >> semilogx(pfa, pd);
%
% This function also computes the Area Under the roc Curve (AUC) on a
% linear scale, where maximum value is 1.0 and a random detector's
% performance would be 0.5 (so if your detector has an AUC metric < 0.5,
% you're better off either doing the opposite of what it tells you or just
% flipping a coin). The routine also computes the AUClg metric, which is
% also an area under the ROC curve, but it emphasizes the low-Pfa portion
% of the ROC space - and is no longer capped to a maximum value of 1.0.
% Rather, AUClg's maximum value is the positive exponent of the smallest
% Pfa... if smallest Pfa measurable was 10^(-5), then AUClg's maximum value
% would be +5. The final output is SCORE, which gives the detector score
% associated with each (Pfa,Pd) point of the ROC curve.
%
% INPUTS:
%   tsTgt  : Vector : Detector test statistic evaluated for target samples
%   tsFa   : Vector : Detector test statistic evaluated for false alarm
%            samples (noise or clutter only)
%
% OUTPUTS:
%   Pfa    : Vector : False alarm probability
%   Pd     : Vector : Detection probability
%   AUC    : Scalar : Area Under the Curve
%   AUClg  : Scalar : Area Under the Curve, but with log-scale on x-axis
%   score  : Vector : Detector score for each ROC curve point
%
% REFERENCE:
%   Fawcett, T. "An Introduction to ROC Analysis," Pattern Recognition
%     Letters, Vol. 27, 2006, pp 861-874. In particular, see Algorithm 1
%     and Algorithm 2 listed in the paper.
%
% MODIFICATION LOG:
%  CMH 10/29/09 Original creation.
%  CMH 04/05/10 Fixed subtle bug having to do with indexing output Pfa & Pd
%                arrays. Added robust check for equality of doubles.
%  CMH 04/08/10 Added computation of score array for each ROC curve point.

[mP,nP] = size(tsTgt); if (mP > nP) tsTgt = tsTgt.'; end  % want row vectors
[mN,nN] = size(tsFa);  if (mN > nN) tsFa  = tsFa.';  end  % want row vectors

P = length(tsTgt);
N = length(tsFa);

L = [ones(1,P) zeros(1,N)];  % L(i) = 1 => P; L(i) = 0 => N;
f = [tsTgt tsFa];
[fsorted,ixSort] = sort(f,2,'descend');
Lsorted = L(ixSort);

FP = 0;     TP = 0;
FPprev = 0; TPprev = 0;
AUC = 0;    AUClg = 0;
fprev = -1e50;
jj = 1;
for ii = 1:(P+N),
  if ( abs(fsorted(ii)-fprev) > eps )
    score(jj) = fsorted(ii);
    Pfa(jj) = FP / N;
    Pd(jj)  = TP / P;
    AUC     = AUC   + TrapezoidArea(FP,FPprev,TP,TPprev);
    if (FP > 1)
      AUClg = AUClg + TrapezoidAreaLogPfa(FP/N,FPprev/N,TP/P,TPprev/P);
    end
    fprev = fsorted(ii);
    FPprev = FP;
    TPprev = TP;
    jj = jj+ 1;
  end
  
  if (Lsorted(ii) == 1)
    TP = TP + 1;
  else
    FP = FP + 1;
  end
end

Pfa(jj)   = FP / N;
Pd(jj)    = TP / P;
score(jj) = fsorted(end) - min(abs(diff(fsorted)));
AUC = AUC + TrapezoidArea(N,FPprev,P,TPprev);
AUC = AUC / (P*N);

AUClg = AUClg + TrapezoidAreaLogPfa(N/N,FPprev/N,P/P,TPprev/P);

end  % function ROC_Curve()


function Area = TrapezoidArea(X1,X2,Y1,Y2)

Base = abs(X1-X2);
AvgHeight = (Y1+Y2)/2;

Area = Base * AvgHeight;

end

function Area = TrapezoidAreaLogPfa(X1,X2,Y1,Y2)

Base = abs(log10(X1)-log10(X2));
AvgHeight = (Y1+Y2)/2;

Area = Base * AvgHeight;

end