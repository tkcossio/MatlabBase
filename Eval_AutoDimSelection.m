% SCRIPT : Eval_AutoDimSelection.m
% PURPOSE: To develop and test an implementation of Zhu and Ghodsi's
% automated approach to finding the knee of a scree plot (plot of
% eigenvalues in descending order). This is a way of doing model selection,
% which Carey Priebe indicates is about the best approach currently in the
% literature for this problem.
%
% REFERENCE:
%  Zhu, M. and A. Ghodsi. "Automatic Dimensionality Selection from the
%    Scree Plot via the use of Profile Likelihood," Computational
%    Statistics & Data Analysis, vol. 51, pgs 918-930, 2006.
%
% MODIFICATION LOG:
%  CMH 02/04/11 Original creation.

% Set up experiment from section 3.1 of reference:
p       = 100;                   % dimension of data
Ntrials = 100;
d_l = 45*rand(Ntrials,50);       % 50 eigenvalues drawn from U[0,45]
d_h = 55 + 45*rand(Ntrials,50);  % 50 eigenvalues drawn from U[55,100]
% d_l = 49*rand(Ntrials,80);       % 80 eigenvalues drawn from U[0,49]
% d_h = 51 + 49*rand(Ntrials,20);  % 20 eigenvalues drawn from U[51,100]
d   = [d_l d_h];                 % total 100 eigenvalues

% Implement automated method of Zhu and Ghodsi:
for ii = 1:Ntrials,
  [qHat(ii),logL(ii,:)] = AutoDimSelect_Zhu(d(ii,:));
end

% Generate representative displays:
figure(1); clf;
[scree,ix] = sort(d(1,:),2,'descend');
plot((1:p),scree,'b-');
xlabel('Dimension');
title('Scree Plot');

figure(2); clf;
plot((1:p),logL(2,:),'r-');
xlabel('Dimension');
title('Profile Log-Likelihood')

med_qHat  = median(qHat);
mad_qHat = mad(qHat,1);
disp(['Estimated dimension: ' num2str(round(100*med_qHat)/100)]);
disp(['Median absolute deviation: (' num2str(round(100*mad_qHat)/100) ')']);
