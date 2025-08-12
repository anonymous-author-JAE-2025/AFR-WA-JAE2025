function [ x ] = logit(x)
x(x==0) = x(x==0)+1e-10;
x(x==1) = x(x==1)-1e-10;
x = log(x./(1-x));
end

