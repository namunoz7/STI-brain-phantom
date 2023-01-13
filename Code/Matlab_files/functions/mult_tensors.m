function tmp1 = mult_tensors(N,At,A)
% Multiply 2 tensors At and A, depending on the size of the tensor. It
% requires the size of the matrix of the images.
% In case the tensors are a matrix, ie, length(size(A)) == 5
if (length(size(A)) == 5)
    tmp1 = zeros([N,size(At,4),size(A,5)]);
    
    for n = 1:size(At,4)
        for m = 1:size(A,5)
            for k = 1:size(At,5)
                tmp1(:,:,:,n,m) = tmp1(:,:,:,n,m) + At(:,:,:,n,k).*A(:,:,:,k,m);
            end
        end
    end
    
% If it is a row of images    
else
    tmp1 = zeros([N,size(At,4)]);
    for n = 1:size(At,4)
        for k = 1:size(At,5)
            tmp1(:,:,:,n) = tmp1(:,:,:,n) + At(:,:,:,n,k).*A(:,:,:,k);
        end
    end
end
end