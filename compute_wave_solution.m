function [F,Fw,V,lambda] = compute_wave_solution(Delta,c, f,df, t)

nt = length(t);


[V,lambda] = svd(-diag(c)*Delta*diag(c));
lambda = diag(lambda);

omega = sqrt(abs(lambda));
fw = V'*(f./c);
dfw = V'*(df./c);
Fw = repmat(fw,[1 nt]) .* cos( omega*t ) + repmat(dfw./omega,[1 nt]) .* sin( omega*t );
cT = repmat(c, [1 nt]);
F = cT .*(V*Fw);