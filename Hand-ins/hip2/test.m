

syms omega

H = 1 / (1 - omega^2 + sqrt(2)*1i*omega);

vpasolve(H*conj(H)==0.5, omega)

vpasolve( abs(H)==1/sqrt(2), omega)
double(solve( abs(H)==1/sqrt(2), omega))


db2mag(-3)

a = [1 2^.5 1];
b = [1];
w = logspace(-1,2,100);
bode(b,a,w)