%program to do interval calculations
function out=interval_stats(U,F,R)

%R=set of stages that count as "reproductive". These are made into new absorbing states

[s,s]=size(U);
Is=eye(s);

%create new transient transition matrix Uprime and mortality matrix Mprime
Uprime=U;
%absorbing state 1 = death before reaching repoductive state
Mprime(1,:)=1-sum(U);
%absorbing state 2 = reproduce before death
%all transitions into R changed to absorption
Mprime(2,:)=sum(Uprime(R,:),1);
%all transitions into R set to 0
Uprime(R,:)=0;

Nprime=inv(eye(s)-Uprime);
Bprime=Mprime*Nprime;
D=diag(Bprime(2,:));


%Fundamental matrix conditional on reaching R before death
Ncond=D*Nprime*inv(D);
%mean time to absorption by reaching R
eta1cond=(sum(Ncond))';
eta2cond=(eta1cond'*(2*Ncond -Is))';
%variance in time to reaching R
varetacond=eta2cond-(eta1cond.*eta1cond);

out.Uprime=Uprime;
out.Mprime=Mprime;
out.Bprime=Bprime;
out.Ncond=Ncond;
out.eta1cond=eta1cond;
out.eta2cond=eta2cond;
out.varetacond=varetacond;
