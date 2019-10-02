%SVM test

rng default

a = readtable('titanic.csv');

surv = table2array(a(:,1));
m1 = table2array(a(:,2));
m2 = table2array(a(:,5:8));
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

data = [m1, m2, genNum];
x = normalize(data);

%high Dim data
% SVMmod = fitcsvm(x, surv);
% 
% CVSVMModel = crossval(SVMmod);
% classLoss = kfoldLoss(CVSVMModel)

%Plain PCA full data
% [V,S,U] = svd(x);
% e1 = U(:,1);
% e2 = U(:,2);
% 
% adj1 = data*e1;
% adj2 = data*e2;
% xp = [adj1, adj2];

%PCA Restricted data
% m2 = table2array(a(:,5:6));
% data = [m2, genNum];
% x = normalize(data);
% [V,S,U] = svd(x);
% e1 = U(:,1);
% e2 = U(:,2);
% 
% adj1 = data*e1;
% adj2 = data*e2;
% xp = [adj1, adj2];

%Gaussian PCA, AGS
% N = 887;
% 
% 
% sigma=.9;
% for i=1:N
%     for j=1:N
%         K(i,j)=exp(-norm(x(i,:)-x(j,:))^2/(2*sigma^2));
%     end
% end
% E=ones(N,N); 
% K_hat=(eye(N,N)-E/N)*K*(eye(N,N)-E/N);
% 
% [V,S]=eig(K_hat); S=diag(S); 
% S=real(S); S=max(S,0); V=real(V);
% S=sqrt(S); [S,I]=sort(S,'descend'); V=V(:,I);
% adj1=S(1)*V(:,1); % first coordinates of the data points
% adj2=S(2)*V(:,2); % second coordinates of the data points

%Poly TESt AGS
% N = 887;
% % PCA with polynomial basis function kernel
% 
% for i=1:N
%     for j=1:N
%         K(i,j)=(dot(x(i,:),x(j,:))^3)/((norm(x(i,:))^3)*(norm(x(j,:))^3));
%     end
% end
% E=ones(N,N); 
% K_hat=(eye(N,N)-E/N)*K*(eye(N,N)-E/N);
% 
% [V,S]=eig(K_hat); S=diag(S); 
% S=real(S); S=max(S,0); V=real(V);
% S=sqrt(S); [S,I]=sort(S,'descend'); V=V(:,I);
% adj1=S(1)*V(:,1); % first coordinates of the data points
% adj2=S(2)*V(:,2);


%Combo KPCA AGS
N = 887;

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

sigma=.9;
for i=1:N
    for j=1:N
        K(i,j)=exp(-norm(x(i,:)-x(j,:))^2/(2*sigma^2));
    end
end
E=ones(N,N); 
K_hat=(eye(N,N)-E/N)*K*(eye(N,N)-E/N);

[V,S]=eig(K_hat); S=diag(S); 
S=real(S); S=max(S,0); V=real(V);
S=sqrt(S); [S,I]=sort(S,'descend'); V=V(:,I);
adj2=S(1)*V(:,1); % first coordinates of the data points



xp = [adj1, adj2];

tic
SVMmod = fitcsvm(xp, surv);
toc

CVSVMModel = crossval(SVMmod);
classLoss = kfoldLoss(CVSVMModel)