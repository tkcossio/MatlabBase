function [pHat,varargout] = AutoDimSelect_Zhu(X)
% CALL:  [pHat] = AutoDimSelect_Zhu(X);
%        [pHat,pLogL] = AutoDimSelect_Zhu(X);
%
% PURPOSE: To apply the automated method of Zhu and Ghodsi 2006 for finding
% the knee of a scree plot, i.e. to determine the number of important
% eigenvalues/dimensions to retain in principal components. The knee is
% estimated via a profile likelihood estimation technique.
%
% INPUTS:
%   X        : Vector : A vector of eigenvalues/singular values, does not
%                        have to be sorted
%
% OUTPUTS:
%   pHat     : Scalar : Estimated number of important principal components,
%                        or estimated dimension.
%
%   Optional:
%   pLogL    : Vector : Estimated profile log-likelihood from which pHat
%                        was estimated.
%
% REFERENCE:
%  Zhu, M. and A. Ghodsi. "Automatic Dimensionality Selection from the
%    Scree Plot via the use of Profile Likelihood," Computational
%    Statistics & Data Analysis, vol. 51, pgs 918-930, 2006.
%
% MODIFICATION LOG:
%  CMH 02/04/11: Original creation.

% Make input vector a row vector:
[m,n] = size(X);
if (m > n) X = X.'; end
p = length(X);  % number of dimensions in input data

% Sort the eigenvalues in order of importance:
xS = sort(X,2,'descend');

logL = NaN*ones(1,p);
for q = 1:p,
  % Partition the p dimensions as [1:q] and [q+1:p]:
  S1 = xS(1:q);
  S2 = xS(q+1:p);

  % Calculate pooled variance estimate:
  s2_1 = std(S1,0,2);
  s2_2 = std(S2,0,2);
  sigSq_hat = ((q-1)*s2_1 + (p-q-1)*s2_2) / (p-2);

  % Calculate means:
  mu_1 = mean(S1,2);
  mu_2 = mean(S2,2);

  % Compute profile likelihood, assuming that eigenvalues are Gaussian
  % distributed with different means and common variance:
  logL(q) = sum(log(normpdf(S1,mu_1,sigSq_hat))) + ...
            sum(log(normpdf(S2,mu_2,sigSq_hat)));
end

% Estimate dimension via maximum of profile log-likelihood:
[pkLogL,pHat] = max(logL);

% Send estimated profile log-likelihood to user, if requested:
if (nargout > 1)
  varargout{1} = logL;
end
