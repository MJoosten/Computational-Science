clear all;
close all;
format compact;
tic;

%variables
N=1000;
K=2000;
l=0.311;%nm
ksi=50;%nm
L=l*K;
sigma = sqrt(l/ksi);
tangent = zeros(2,K,N);
tangent(1,1,:)=1;
shape = zeros(2,K,N);
shape(1,1,:)=1;
unittangent = zeros(2,K,N);
Rsq = zeros(N,1);
angles = zeros(K,N);

for nn=1:N
    for kk=2:K
        psi = sigma.*randn(1);
        angles(kk,nn) = psi;
        tangent(:,kk,nn) = [cos(psi) -sin(psi);sin(psi) cos(psi)]*tangent(:,kk-1,nn);
        unittangent(:,kk,nn) = tangent(:,kk,nn)./norm(tangent(:,kk,nn));
        shape(:,kk,nn) = shape(:,kk-1,nn)+l*unittangent(:,kk,nn);
    end
    Rsq(nn) = norm((shape(:,end,nn)-shape(:,1,nn))).^2;
    avRsq(nn) = sum(Rsq(1:nn))./nn;
    stdRsq(nn) = (std(Rsq(1:nn)))./sqrt(nn);
end
time = toc;

RsqEx = 4*ksi*L-(8*ksi^2)*(1-exp(-L/2*ksi));
%% plots & prints
figure(1)
subplot(2,2,1)
for nn=1:N/10
    colors = (1/N)*randperm(N,3);
    quiver(shape(1,:,nn),shape(2,:,nn),unittangent(1,:,nn),unittangent(2,:,nn),'color',colors);hold on;
    title('plot of DNA in random configuration');xlabel('x [nm]');ylabel('y [nm]');
end
subplot(2,2,2)
quiver(shape(1,:,1),shape(2,:,1),unittangent(1,:,1),unittangent(2,:,1))
title('sinle strand plotted');xlabel('x [nm]');ylabel('y [nm]');
subplot(2,2,3)
histogram(angles(:,1));title('distribution of bend angles');xlabel('angle [radians]');ylabel('count')
subplot(2,2,4)
plot(1:N,avRsq,1:N,stdRsq,1:N,RsqEx*ones(N,1));title('convergence of R-distance with number of simulated polymers');xlabel('number of poymers');ylabel('R-distance [nm^2]')
legend('average R-distance','standard deviation of R-distance','expected R-distance')
sprintf('average R-distance = %0.3fnm^2 \nstandard deviation R-distance = %0.3fnm^2',avRsq(end),stdRsq(end))
sprintf('expected R-distance = %0.3f',RsqEx)
sprintf('computing time %0.3fseconds',time)