%program to calculate statistics of reproductive events
function out=repro_stats(U,F);

[s,s]=size(U);
Is=eye(s);

%successful reproductive events
R=2;

%moments of occupancy times
%first moments
N1 = inv(Is-U);
%second moments
N2 = (2*diag(diag(N1))-Is)*N1;
%variance of occupancy
varN=N2-(N1.*N1);
%sd of occupancy
sdN=sqrt(varN);
%cv of occupancy
cvN=sdN./N1;

%net reproductive rate
R0=max(eig(F*N1));

out.N1=N1;
out.N2=N2;
out.varN=varN;
out.sdN=sdN;
out.cvN=cvN;
out.R0=R0;



