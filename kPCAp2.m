%Code to generate figure 6
a = readtable('titanic.csv');

surv = table2array(a(:,1));
m1 = table2array(a(:,2));
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
x = normalize(data);
N = 887;
% PCA with polynomial basis function kernel

for i=1:N
    for j=1:N
        K(i,j)=(dot(x(i,:),x(j,:))^3)/((norm(x(i,:))^3)*(norm(x(j,:))^3));
    end
end
E=ones(N,N); 
K_hat=(eye(N,N)-E/N)*K*(eye(N,N)-E/N);

[V,S]=eig(K_hat); S=diag(S); 
S=real(S); S=max(S,0); V=real(V);
S=sqrt(S); [S,I]=sort(S,'descend'); V=V(:,I);
adj1=S(1)*V(:,1); % first coordinates of the data points
adj2=S(2)*V(:,2); % second coordinates of the data points

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
title('Polynomial KPCA, only Age, Gender, Siblings');
shg;

