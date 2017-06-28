%% Computational Science Final Project: Worm-Like Chain
% Task 4
% Authors: Maarten Joosten & Nemo Andrea
% IDs: xxxxxxx & 4473035
% Date of Creation: 26-06-2017
% github: https://github.com/MJoosten/Computational-Science

clear all;
close all;
format compact;

%% Start

%----------------------TODO: VECTORISE, BABY ------------------------------
%parameters
P=8; %number of configurations (configuration = amount segments of chain) 
P_range=[10,5000]; %range of segment numbers (default: 100,5000)
N=100; %Iterations of Polymer/chain (DNA) generation (default:100)
K=round(linspace(P_range(1),P_range(2),P)); % Number of segments of chain
                                     %(base pairs) (default:2000)
length_link=0.311;%[nm] Length of each chain link(base pair)(default:0.311)
length_persist=50; %[nm] persistence length (default:50)
length_chain=K*length_link; %[nm] Total length of chain (DNA)
t_initial=[0;0;1]; %initial orientation of t vector (unit length);
                   %(default: 0,0,1 (z axis))
Nbins = 30;
Npoints = 300;
sigma = sqrt(length_persist.*length_chain./2);
r = 3;

%Preallocation - Outside Loop
comp_time=zeros(N,1); %this array will hold the computational time for
                      %for each Iteration (N iterations)
distances=zeros(N,P); %will hold the squared end-to-end distances
    Xend = zeros(P,N);
    Yend = zeros(P,N);
    bins = zeros(P,Nbins);
    Px = zeros(P,Npoints);
    points = zeros(P,Npoints);
    devX = zeros(P,1);
    devY = zeros(P,1);
    
%opening statement (for console iterpretability)
fprintf(['\n>>>[task 3] Starting Computation WLC 3D with %u'...
        ' configurations each with %u iterations and number of segments '...
        'between %u and %u'],P,N,min(K),max(K))
for pp=1:P
    K_local=K(pp);
    N_local=N;
    %Preallocation - Inside Loop
    location=zeros(3,K_local,N_local); %will hold the location for each polymer link (3D)
    tangents=ones(3,K_local,N_local);% holds the angles %TODO: do this more efficiently
    tangents(1,:,:)=tangents(1,:,:)*t_initial(1); %setting initial tangent
    tangents(2,:,:)=tangents(2,:,:)*t_initial(2); %setting initial tangent
    tangents(3,:,:)=tangents(3,:,:)*t_initial(3); %setting initial tangent

    % generate random bend angles
    % Gaussian Distribution with mu=0;var=length_link/length_persistence
    rand_angles=sqrt(length_link/length_persist)*randn(2,K_local,N_local);

    % Computation ------------------------------------------------------------- 

    fprintf('\nComputing WLC 3D Distance for K=%u links, for N=%u iterations',K_local,N_local)
    for ii=1:N_local %loop over N iterations(generate N independent runs)
        tic %start a clock for each run 
        for jj=1:K_local-1 %compute K segments %FIX                      
            %find alpha and beta of PREVIOUS iteration
            alpha_t=acos(tangents(3,jj,ii)); %arccos(t_z)       
            beta_t=atan2(tangents(2,jj,ii),tangents(1,jj,ii));%arctan(t_y/t_x)
            ortho_1=[cos(alpha_t)*cos(beta_t);cos(alpha_t)*sin(beta_t);-sin(alpha_t)];
            ortho_2=[-sin(beta_t);cos(beta_t);0];

            %calculate coefficients
            %TODO: i feel like we can do some clever rewriting to reduce
            %computational cost here (trigonometry???)
            norm_factor=sqrt(1-(sin(rand_angles(1,jj,ii))*sin(rand_angles(2,jj,ii)))^2);
            c_t=(cos(rand_angles(1,jj,ii))*cos(rand_angles(2,jj,ii)))/norm_factor;
            c_1=(sin(rand_angles(1,jj,ii))*cos(rand_angles(2,jj,ii)))/norm_factor;
            c_2=(cos(rand_angles(1,jj,ii))*sin(rand_angles(2,jj,ii)))/norm_factor; 

            %calculate the new tangent vector (3D)
            tangents(:,jj+1,ii)=c_t*tangents(:,jj,ii)+c_1*ortho_1+c_2*ortho_2;    
        end
        
        %update Locations (fast method)
        location(:,:,ii)=cumsum(tangents(:,:,ii)*length_link,2); 
    end
    
    Xend(pp,:) = location(1,end,:);
    Yend(pp,:) = location(2,end,:);
    bins(pp,:) = linspace(-r*sigma(pp),r*sigma(pp),Nbins);
    points(pp,:) = linspace(-r*sigma(pp),r*sigma(pp),Npoints);
    %theoretical distribution
    Px(pp,:) = 1./sqrt(pi*length_persist*length_chain(pp))*exp((-points(pp,:).^2)./(length_persist*length_chain(pp)));
end



%signaling computation is finished
fprintf('\n>%u Configurations each with %u iterations completed, Computation finished',P,N)
%% plots
   figure
for pp=1:P
    subplot(2,P,pp)
    Hx = hist(Xend(pp,:),bins(pp,:));
    bar(bins(pp,:),Hx);
    title(sprintf('histogram of end position in x over %i iterations and %i segments',N,K(pp)));xlabel('position [nm]');ylabel('count')
    subplot(2,P,pp+P)
    Hy = hist(Yend(pp,:),bins(pp,:));hold on;
    bar(bins(pp,:),Hy);
    title(sprintf('histogram of end position in y over %i iterations and %i segments',N,K(pp)));xlabel('position [nm]');ylabel('count') 

%     sprintf('standard deviation in x for K=%i: %0.4f \nmean in x for K=%i: %0.4f \npredicted deviation: %0.4f',K(pp),std(Xend(pp,:)),K(pp),mean(Xend(pp,:)),sigma(pp))
%     sprintf('standard deviation in y for K=%i: %0.4f \nmean in y for K=%i: %0.4f \npredicted deviation: %0.4f',K(pp),std(Yend(pp,:)),K(pp),mean(Yend(pp,:)),sigma(pp))

    %calculate the standard deviations of respective coordinates for each value of K
    devX(pp)=std(Xend(pp,:));
    devY(pp)=std(Yend(pp,:));
end


figure
subplot(1,2,1)
bar(bins(end,:),Hx);title(sprintf('zoom into K=%i segments (x)',K(end)))
subplot(1,2,2)
bar(bins(end,:),Hy);title(sprintf('zoom into K=%i segments (y)',K(end)))

figure
subplot(1,2,1)
plot(K,devX,K,sigma);title('standard deviation in x as a function of segment length');grid on;
legend('simulated deviation','predicted deviation')
subplot(1,2,2)
plot(K,devY,K,sigma);title('standard deviation in y as a function of segment length');grid on;
legend('simulated deviation','predicted deviation')