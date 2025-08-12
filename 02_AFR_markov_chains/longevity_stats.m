%longevity stats
function out=longevity_stats(U)
%program to calculate statistics of longevity from U matrix

[s,s]=size(U);
Is=eye(s);

%survivorship 
n0=Is(:,1);
tlimit=100;
for t=0:tlimit
    n(:,t+1)=(U^t)*n0;
end
survivorship=sum(n)';

%fundamental matrix
N1=inv(Is-U);

%statistics of longevity
eta1=sum(N1)';
eta2=(eta1'*(2*N1-Is))'; % second moments
vareta=eta2-(eta1.*eta1);
sdeta=sqrt(vareta);
cveta=sdeta./eta1;

out.U=U;
out.survivorship=survivorship;
out.N1=N1;
out.eta1=eta1;
out.eta2=eta2;
out.vareta=vareta;
out.sdeta=sdeta;
out.cveta=cveta;

