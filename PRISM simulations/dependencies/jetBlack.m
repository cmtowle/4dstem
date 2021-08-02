function [cmap] = jetBlack(Nc)

if nargin < 1
    Nc = 1024;
end

cmap = jet(Nc);
b = round(Nc/8);
c = linspace(0,1,b)';
cmap(1:b,:) = cmap(1:b,:) .* c;

end

