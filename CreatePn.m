function pn = CreatePn(N,M)


mseq = round(rand(N));
pn = zeros(N,1);
if (nargin < 1)&&(nargin > 2)
    error('Wrong number of input arguments');
    
elseif (nargin == 1)
    pn = round(rand(1,N))*2-1;
    %pn = (mseq(1:N)*2-1)';
    
elseif (nargin == 2)
    rate = round(N/M);
    an = round(rand(1,M))*2-1;
    %an = mseq(1:M)*2-1;
    for i = 1:M
       pn((i-1)*rate+1:i*rate) = an(i);  
    end
end