function [alpha1, alpha2, alpha3] = dfapiece(d,p)

[m,~] = size(d);
tempY1 = [];
tempX1 = [];
tempY2 = [];
tempX2 = [];
tempY3 = [];
tempX3 = [];

for i = 1:m
    if (d(i) >= 10) && (d(i) <= 30) % 10<=n<=30
        tempY1(end+1) = p(i);
        tempX1(end+1) = d(i);
    end
    if (d(i) >= 30) && (d(i) <= 300) % 30<=n<=300
        alpha1 = polyfit(log10(tempX1),log10(tempY1),1);
        alpha1 = alpha1(1);
        tempY2(end+1) = p(i);
        tempX2(end+1) = d(i);
    end
    if (d(i) >= 300) % 300<=n
        alpha2 = polyfit(log10(tempX2),log10(tempY2),1);
        alpha2 = alpha2(1);
        tempY3(end+1) = p(i);
        tempX3(end+1) = d(i);
    end
end

alpha3 = polyfit(log10(tempX3),log10(tempY3),1);
alpha3 = alpha3(1);