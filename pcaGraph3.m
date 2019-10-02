%Code to generate figure 3
a = readtable('titanic.csv');

surv = table2array(a(:,1));
m2 = table2array(a(:,5:6));

genChar = table2array(a(:,4));

genNum = [];
for i = 1:887
    temp = 1;
    if strcmp(genChar(i), 'male')
        genNum(i) = 1;
    else
        genNum(i) = 0;
    end
end
genNum = genNum';

data = [m2, genNum];
data = normalize(data);
[V,S] = svd(data);
e1 = V(:,1);
e2 = V(:,2);

adj1 = data*e1;
adj2 = data*e2;
figure(1);
subplot(111);
hold on;

for i = 1:887
    
    if surv(i) == 1
        tmp = '.b';
    else
        tmp = '.r';
    end
    plot(adj1(i), adj2(i), tmp)
end
hold off;
shg;

