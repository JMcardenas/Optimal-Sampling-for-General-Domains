function f = F_Genz(y)
%-------------------------------------------------------------------------%
% Description: function for compute the Genz product peak
% Programer: Juan Manuel Cardenas
% Date: May 10-2019 / Last modification : 
%-------------------------------------------------------------------------%
y = y';
[d,m] = size(y);
ind = [1:d]';
aux = ind + 1;
v_1 = (y + (( (-1).^aux )./aux ) ).^2 ; 
v_2 = (d/4) + v_1 ; 
v_3 = (d/4)./(v_2);
f   = prod(v_3);
end
