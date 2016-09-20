function y = compute_laplacian(x, c)


n = length(c);
c = repmat(c, [1, size(x,2)]);
w = 4*pi^2 * [0:n/2 -n/2+1:-1]'.^2;
w = repmat(w, [1, size(x,2)]);

y = real( ifft(fft(c.*x,[],1).*w,[],1).*c );