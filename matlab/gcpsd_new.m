function [S, w] = gcpsd_new(input)
% -------------------------------------------------------------------------
%   GCPSD computes the Generalized Cross Power Spectral Density
% -------------------------------------------------------------------------

% Check for correct dimensions
if size(input, 1) < size(input, 2)
    input = input';
end

[len,NoC] = size(input);
fft_len = max(round(len/4 + 1), 2^7 + 1);
Px1 = pwelch(input,len,[],2*(fft_len-1) + 1);
Px2 = pwelch(input,round(len/2),[],2*(fft_len-1) + 1);
Px3 = pwelch(input,round(len/4),[],2*(fft_len-1) + 1);
Px4 = pwelch(input,round(len/6),[],2*(fft_len-1) + 1);
% Px = Px1.*Px2.*Px3.*Px4;
Px = Px1 + Px2 + Px3 + Px4;

for i = 1:size(Px,2)
    sigma(:,:,i) = Px(:,i).*Px;
end

sigma(:,eye(NoC)==1) = 1;
sigma = permute(sigma, [2 3 1]);
for i=1:fft_len
    lmax(i) = eigs(sigma(:, :, i), 1);
end
S = ((lmax - 1) ./ (NoC - 1)).^2;
w = (0:pi/(fft_len-1):pi).';
if size(S,1) < size(S,2)
    S = S';
end

% % Initialize constants 
% [len, NoC] = size(input);
% fft_len = max(round(len/4 + 1), 2^7 + 1);
% lmax  = zeros(fft_len, 1);
% 
% % Compute pairwise cross power spectral densities for lower triangular
% %---------------------Charalampos (raplaced fors)--------------------------
% sigma = abs(cpsd(input, input, window, [], 2*(fft_len-1) + 1,'mimo'));
% sigma(:,eye(NoC)==1) = 1;
% %--------------------------------------------------------------------------
% 
% % Permute matrix for computational purposes
% sigma = permute(sigma, [2 3 1]);
% 
% %---------------------Charalampos--------------------------
% for i=1:fft_len
%     lmax(i) = eigs(sigma(:, :, i), 1);
% end
% %--------------------------------------------------------------------------
% 
% % GCPSD formula
% S = ((lmax - 1) ./ (NoC - 1)).^2;
% w = (0:pi/(fft_len-1):pi).';

end