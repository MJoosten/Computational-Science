clear all;
close all;
format compact;
tic;

%variables
N=100;
K=2000;
l=0.311;%nm
ksi=50;%nm
sigma = sqrt(l/ksi);
tangent = zeros(2,K,N);
tangent(1,1,:)=1;
shape = zeros(2,K,N);
shape(1,1,:)=1;
unittangent = zeros(2,K,N);
Rsq = zeros(N,1);

for nn=1:N
    for kk=2:K
        psi = sigma.*randn(1);
        tangent(:,kk,nn) = [cos(psi) -sin(psi);sin(psi) cos(psi)]*tangent(:,kk-1,nn);
        unittangent(:,kk,nn) = tangent(:,kk,nn)./norm(tangent(:,kk,nn));
        shape(:,kk,nn) = shape(:,kk-1,nn)+l*unittangent(:,kk,nn);
    end
    Rsq(nn) = norm((shape(:,end,nn)-shape(:,1,nn))).^2;
end
avRsq = sum(Rsq)./N;
stdRsq = (std(Rsq))./sqrt(N);
time = toc;

%% plots & prints
figure(1)
for nn=1:N
    colors = (1/N)*randperm(N,3);
    quiver(shape(1,:,nn),shape(2,:,nn),unittangent(1,:,nn),unittangent(2,:,nn),'color',colors);hold on;
end
figure(2)
quiver(shape(1,:,1),shape(2,:,1),unittangent(1,:,1),unittangent(2,:,1))

sprintf('average R-distance = %0.3fnm^2 \nstandard deviation R-distance = %0.3fnm^2',avRsq,stdRsq)
sprintf('computing time %0.3fseconds',time)