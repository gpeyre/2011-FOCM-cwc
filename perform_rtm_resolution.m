function R = perform_rtm_resolution(Ffwd, Fback, c0)

% R = perform_rtm_resolution(Ffwd, Fback, c0);

nt = size(Fback,2);
n = size(Fback,1);

% compute derivative in time
sel1 = [1 1:nt-1];
sel2 = [2:nt nt];
Ffwd = ( 2*Ffwd - Ffwd(:,sel1) - Ffwd(:,sel2) )*nt^2;

% RTM solution
R = - 2./c0 .* sum( Ffwd.*Fback, 2 )/nt;
selz = [1:round(.3*n), round(.9*n):n];
R(selz) = 0;