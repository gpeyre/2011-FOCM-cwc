function df = compute_spectral_derivative(f, use_spectral)

if nargin<2
    use_spectral=1;
end

n = size(f,1);

if use_spectral
    x = [0:n/2 -n/2+1:-1]';
    p = size(f,2);
    x = x*ones(1,p);
    df = ifft( fft(f) .* (2*pi*1i*x) );
    
%    for i=1:size(f,2)
%        progressbar(i,size(f,2));
%        df(:,i) = ifft( fft(f(:,i)) .* (2*pi*1i*x) );
%    end
    df = real(df);
else
    df = (f([2:end,1],:)-f([end 1:end-1],:))*n/2;
end