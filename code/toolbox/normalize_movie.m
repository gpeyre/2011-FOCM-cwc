function Adisp = normalize_movie(A)

% normalize_movie - normalie movie file for display

Adisp = A;
for i=1:size(A,3)
    Adisp(1,1,i) = -max(max(abs(A(:,:,i)))); 
    Adisp(2,1,i) = max(max(abs(A(:,:,i)))); 
end