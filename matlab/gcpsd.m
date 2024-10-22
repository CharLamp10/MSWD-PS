function [S, w] = gcpsd(input, window)
% -------------------------------------------------------------------------
%   GCPSD computes the Generalized Cross Power Spectral Density
% -------------------------------------------------------------------------

if nargin<2
    window = [];
end

% Check for correct dimensions
if size(input, 1) < size(input, 2)
    input = input';
end

% Initialize constants 
[len, NoC] = size(input);
fft_len = max(round(len/4 + 1), 2^7 + 1);
sigma = zeros(fft_len, NoC, NoC);
lmax  = zeros(fft_len, 1);

% Compute pairwise cross power spectral densities for lower triangular
% matrix (since it is symmetric)

% for i=1:NoC
%     for j=1:NoC
%         if i == j
%             sigma(:, i, j) = ones(fft_len, 1);
%         elseif i > j
%             sigma(:, i, j) = abs(cpsd(input(:, i), input(:, j), window, [], 2*(fft_len-1) + 1));
%         end
%     end
% end
%---------------------Charalampos (raplaced fors)--------------------------
sigma = abs(cpsd(input, input, window, [], 2*(fft_len-1) + 1,'mimo'));
sigma(:,eye(NoC)==1) = 1;
%--------------------------------------------------------------------------
% Permute matrix for computational purposes
sigma = permute(sigma, [2 3 1]);

% Fill symmetric matrix and calculate max eigenvalue for every omega
%---------------------Charalampos--------------------------
for i=1:fft_len
    %sigma(:, :, i) = sigma(:, :, i) + tril(sigma(:, :, i), -1).'; No need
    %after the change with the fors
    lmax(i) = eigs(sigma(:, :, i), 1);
end
%--------------------------------------------------------------------------

% GCPSD formula
S = ((lmax - 1) ./ (NoC - 1)).^2;
w = (0:pi/(fft_len-1):pi).';

end
