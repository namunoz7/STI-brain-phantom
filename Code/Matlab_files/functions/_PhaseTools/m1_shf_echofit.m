function [lfs,res,off] = m1_shf_echofit(ph,mag,TE,inter)
% ECHOFIT Magnitude-weighted least square regression of phase to echo time
%   [LFS,RES] = echofit(PH,MAG,TE,INTER)
%
%   PH:    unwrapped phases from multiple echoes
%   MAG:   corresponding magnitudes 
%   TE:    echo times
%   INTER: intercept of linear fitting, zero or non-zero
%   LFS:   local field shift after fitting
%   RES:   fitting residuals
%   OFF:   constant offset fitting intercept
% 
% Adapted from original code by HF Sun
% Last modified by Julio Acosta-Cabronero

% check ph and mag have same dimensions
if ~size(ph)==size(mag)
    error('ERROR: Phase and magnitude matrices not the same size');
end

if ~exist('inter','var') || isempty(inter)
    inter = 0; % default: zero intercept
end

[np,nv,ns,ne] = size(ph);

ph  = permute(ph, [4 1 2 3]);
mag = permute(mag,[4 1 2 3]);

ph  = reshape(ph,ne,[]);
mag = reshape(mag,ne,[]);

if ~inter % zero intercept
	TE_rep = repmat(TE(:),[1 np*nv*ns]);

	lfs = sum(mag.*ph.*TE_rep,1)./(sum(mag.*TE_rep.*TE_rep)+eps);
	lfs = reshape(lfs,[np nv ns]);

	% calculate fitting residual
	lfs_rep = permute(repmat(lfs(:),[1 ne]),[2 1]);
	res = reshape(sum((ph - lfs_rep.*TE_rep).*mag.*(ph - lfs_rep.*TE_rep),1)./sum(mag,1)*ne,[np nv ns]);% % normalize lfs ("/sum*ne")
	res(isnan(res)) = 0; res(isinf(res)) = 0;
    
else % non-zero intercept
	x = [TE(:), ones(length(TE),1)];
	beta = zeros(2, np*nv*ns);
	res = zeros([np nv ns]);
	matlabpool open
	parfor i = 1:np*nv*ns
		y = ph(:,i);
		w = mag(:,i);
		beta(:,i) = (x'*diag(w)*x)\(x'*diag(w)*y);
		res(i) = (y-x*beta(:,i))'*diag(w)*(y-x*beta(:,i))/sum(w)*ne;
	end
	matlabpool close
	beta(isnan(beta)) = 0;
	beta(isinf(beta)) = 0;
	res(isnan(res)) = 0;
	res(isinf(res)) = 0;

	lfs = reshape(beta(1,:),[np nv ns]);
	off = reshape(beta(2,:),[np nv ns]);
	res = reshape(res,[np nv ns]);

end

% JAC
if length(TE)>1
    lfs = lfs*(TE(2)-TE(1));
else
    lfs = lfs*TE;
end
    
% % create sparse x matrix
% x = sparse(length(TE)*np*nv*ns,2*length(TE)*np*nv*ns);
% x(sub2ind([length(TE)*np*nv*ns,2*length(TE)*np*nv*ns],1:length(TE)*np*nv*ns,1:2:2*length(TE)*np*nv*ns)) = 1;
% x(sub2ind([length(TE)*np*nv*ns,2*length(TE)*np*nv*ns],1:length(TE)*np*nv*ns,2:2:2*length(TE)*np*nv*ns)) = repmat(TE,[np*nv*ns,1]);

% beta = inv(x'*diag(mag(:))*x)*(x'*(mag(:).*ph(:)));
% beta = reshape(beta,2,:);
% lfs = reshape(beta(1,:),np,nv,ns);
% off = reshape(beta(2,:),np,nv,ns);
