function V_filt = filt_eigenvectors(V)
% Filters the DTI eigenvectors using a median filter
V_filt = zeros(size(V));
% median_kernel = [5,5,5; 3,3,3; 3,3,3];
median_kernel = 3*ones(3);
for nn=1:3
    for mm=1:3
        disp([nn, mm])
        V_filt(:, :, :, nn, mm) = medfilt3(V(:, :, :, nn, mm), median_kernel(nn,:));
    end
end
end
