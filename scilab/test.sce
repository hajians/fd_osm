// For reading a sparse matrix

function Asp = readsparse(filename)
// readsparse reads a sparse matrix made by a FORTRAN code and
// gives a sparse matrix in scilab.
  init = fscanfMat(filename);
  dim  = zeros(2,1)
  dim(1) = init(1,1);
  dim(2) = init(1,2);
  
  len  = init(1,3)+1;
  
  Asp = sparse([init(2:len,1),init(2:len,2)], init(2:len,3), dim)
endfunction

stacksize('max');

A1    = readsparse('A1');
A1GAM = readsparse('A1GAM');
AGAM  = readsparse('AGAM');
N     = size(A1,1);

A2    = readsparse('A2');
A2GAM = readsparse('A2GAM');

B1 = (A1GAM'/A1)*A1GAM;
B2 = (A2GAM'/A2)*A2GAM;

AGAMsq    = full(AGAM)^0.5;
invAGAMsq = inv(AGAMsq);

C1 = invAGAMsq * B1 * invAGAMsq;
C2 = invAGAMsq * B2 * invAGAMsq;

D1 = C1 / ( speye(C1) - C1 );
invD1 = speye(D1)/D1;

h = 1.0/sqrt(N);
Gamma = 0.5*(1+h^(0/2));
D2 = - speye(C2) + (2*Gamma-1) / (Gamma*speye(C2)-C2);

//ev = spec(full(B1),full(AGAM));
//ev1 = spec(C2);

//ev2 = spec(D2);
ev3 = spec(D1);
Maxev = max(real(ev3));
Minev = min(real(ev3));

    
printf('%d %f %f',sqrt(N), Maxev, Minev); 