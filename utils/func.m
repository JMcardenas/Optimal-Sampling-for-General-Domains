function y = func(varagin)

z=varagin';
[d,m]=size(z);
y=exp(-sum(z)/(d)); % dim comp function
y=y';

end

