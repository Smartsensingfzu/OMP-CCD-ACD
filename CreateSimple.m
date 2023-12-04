function a = CreateSimple(N,M)
%产生采样矩阵M*N

if (nargin ~= 2)
    error('Wrong number of input arguments');
    
else
    a = zeros(M,N);
    rate = floor(N/M);
    for i = 1:M
        for j = 1:rate
           a(i,(i-1)*rate+j) = 1;
        end
    end
end