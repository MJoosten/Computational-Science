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

%parameters
P=20; %number of configurations (configuration = amount segments of chain) 
P_range=[100,5000]; %range of segment numbers (default: 100,5000)
S = 5; % number of runs (runs = number of iterations)
S_range = [10 1000]; %range of iteration number
defN=100; %Iterations of Polymer/chain (DNA) generation (default:100)
K=round(linspace(P_range(1),P_range(2),P)); % Number of segments of chain
                                     %(base pairs) (default:2000)
defK = 1000; %dfault number of segments
N = round(linspace(S_range(1),S_range(2),S)); % iteration number
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
comp_time=0; %this array will hold the computational time for
                      %for each Iteration (N iterations)
distances=zeros(defN,P); %will hold the squared end-to-end distances
    Xend = zeros(P,defN);
    Yend = zeros(P,defN);
    bins = zeros(P,Nbins);
    Px = zeros(P,Npoints);
    points = zeros(P,Npoints);
    devX = zeros(P,1);
    devY = zeros(P,1);
    
%opening statement (for console iterpretability)
fprintf(['\n>>>[task 4] Starting Computation WLC 3D with %u'...
        ' configurations each with %u iterations and number of segments '...
        'between %u and %u'],P,defN,min(K),max(K))
for pp=1:P
     K_local=K(pp);    
     
    %Preallocation - Inside Loop
    location=zeros(3,K_local,defN); %will hold the location for each polymer link (3D)
    tangents=ones(3,K_local,defN);% holds the angles 
    norm_factor=zeros(defN,1);
    ortho_1=zeros(defN,1);
    ortho_2=zeros(defN,1);
    alpha_t=zeros(defN,1);
    beta_t=zeros(defN,1);
    c_t=zeros(defN,1);
    c_1=zeros(defN,1);
    c_2=zeros(defN,1);
    
    %TODO: do this more efficiently
    tangents(1,:,:)=tangents(1,:,:)*t_initial(1); %setting initial tangent
    tangents(2,:,:)=tangents(2,:,:)*t_initial(2); %setting initial tangent
    tangents(3,:,:)=tangents(3,:,:)*t_initial(3); %setting initial tangent

    
    % generate random bend angles
    % Gaussian Distribution with mu=0;var=length_link/length_persistence
     rand_angles=sqrt(length_link/length_persist)*randn(2,K_local,defN);
     cos_1=reshape(cos(rand_angles(1,:,:)),[K_local,defN]);
     sin_1=reshape(sin(rand_angles(1,:,:)),[K_local,defN]);
     cos_2=reshape(cos(rand_angles(2,:,:)),[K_local,defN]);
     sin_2=reshape(sin(rand_angles(2,:,:)),[K_local,defN]);
     
    % Computation ------------------------------------------------------------- 

    fprintf('\nComputing WLC 3D Distance for K=%u links, for N=%u iterations',K_local,defN)
    tic %start a clock for each run 
           
    
    for jj=1:K_local-1 %compute K segments %FIX                      
        %find alpha and beta of PREVIOUS iteration
        alpha_t=reshape(acos(tangents(3,jj,:)),[1,defN]); %arccos(t_z)       
        beta_t=reshape(atan2(tangents(2,jj,:),tangents(1,jj,:)),[1,defN]);%arctan(t_y/t_x)            
       
        ortho_1=[cos(alpha_t).*cos(beta_t);cos(alpha_t).*sin(beta_t);-sin(alpha_t)];
        ortho_2=[-sin(beta_t);cos(beta_t);zeros(1,defN)];

        %calculate coefficients       
        norm_factor=sqrt(1-(sin_1(jj,:).*sin_2(jj,:)).^2);
        c_t=(cos_1(jj,:).*cos_2(jj,:))./norm_factor;
        c_1=(sin_1(jj,:).*cos_2(jj,:))./norm_factor;
        c_2=(cos_1(jj,:).*sin_2(jj,:))./norm_factor;        

        %calculate the new tangent vector (3D)
        tangents(1,jj+1,:)=c_t.*reshape(tangents(1,jj,:),[1,defN])+c_1.*ortho_1(1,:)+c_2.*ortho_2(1,:);
        tangents(2,jj+1,:)=c_t.*reshape(tangents(2,jj,:),[1,defN])+c_1.*ortho_1(2,:)+c_2.*ortho_2(2,:);
        tangents(3,jj+1,:)=c_t.*reshape(tangents(3,jj,:),[1,defN])+c_1.*ortho_1(3,:)+c_2.*ortho_2(3,:);
    end

    %update Locations (fast method)
    location=cumsum(tangents*length_link,2); 

    %Compute the squared end-to-end distance (works for non-zero starting
    %points too. Alternative method would be norm(vector)^2.      
    distances(:,pp)=sum((location(:,end,:)-location(:,1,:)).^2);
     
    
    comp_time=toc; %clock in computation time for this  WLC set. 
    
    Xend(pp,:) = location(1,end,:);
    Yend(pp,:) = location(2,end,:);
    bins(pp,:) = linspace(-r*sigma(pp),r*sigma(pp),Nbins);
    points(pp,:) = linspace(-r*sigma(pp),r*sigma(pp),Npoints);
    %theoretical distribution
    Px(pp,:) = 1./sqrt(pi*length_persist*length_chain(pp))*exp((-points(pp,:).^2)./(length_persist*length_chain(pp)));
end
fprintf('\n>%u Configurations each with %u iterations completed, Computation 1 finished',P,defN)   
%%
distances2=zeros(defK,S); %will hold the squared end-to-end distances
Xend2 = zeros(S,defK);
Yend2 = zeros(S,defK);
devX2 = zeros(S,1);
devY2 = zeros(S,1);
length_chain2 = defK*length_link;
sigma2 = sqrt(length_persist*length_chain2/2);
bins2 = linspace(-r*sigma2,r*sigma2,Nbins);
points2 = linspace(-r*sigma2,r*sigma2,Npoints);
%theoretical distribution
Px = 1./sqrt(pi*length_persist*length_chain2)*exp((-points2.^2)./(length_persist*length_chain2));

for ss=1:S
    K_local=defK;
    N_local=N(ss);
    %Preallocation - Inside Loop
    location2=zeros(3,K_local,N_local); %will hold the location for each polymer link (3D)
    tangents2=ones(3,K_local,N_local);% holds the angles 
    %TODO: do this more efficiently
    tangents2(1,:,:)=tangents2(1,:,:)*t_initial(1); %setting initial tangent
    tangents2(2,:,:)=tangents2(2,:,:)*t_initial(2); %setting initial tangent
    tangents2(3,:,:)=tangents2(3,:,:)*t_initial(3); %setting initial tangent
   
    % generate random bend angles
    % Gaussian Distribution with mu=0;var=length_link/length_persistence
    rand_angles=sqrt(length_link/length_persist)*randn(2,K_local,N_local);   
    
    % Computation ------------------------------------------------------------- 

    fprintf('\nComputing WLC 3D Distance for K=%u links, for N=%u iterations',K_local,N_local)
    tic %start a clock for each run 
    for ii=1:N_local %loop over N iterations(generate N independent runs)
        for jj=1:K_local-1 %compute K segments %FIX                      
            %find alpha and beta of PREVIOUS iteration
            alpha_t=acos(tangents2(3,jj,ii)); %arccos(t_z)       
            beta_t=atan2(tangents2(2,jj,ii),tangents2(1,jj,ii));%arctan(t_y/t_x)
            ortho_1=[cos(alpha_t)*cos(beta_t);cos(alpha_t)*sin(beta_t);-sin(alpha_t)];
            ortho_2=[-sin(beta_t);cos(beta_t);0];

            %calculate coefficients           
            norm_factor=sqrt(1-(sin(rand_angles(1,jj,ii))*sin(rand_angles(2,jj,ii)))^2);
            c_t=(cos(rand_angles(1,jj,ii))*cos(rand_angles(2,jj,ii)))/norm_factor;
            c_1=(sin(rand_angles(1,jj,ii))*cos(rand_angles(2,jj,ii)))/norm_factor;
            c_2=(cos(rand_angles(1,jj,ii))*sin(rand_angles(2,jj,ii)))/norm_factor; 

            %calculate the new tangent vector (3D)
            tangents2(:,jj+1,ii)=c_t*tangents2(:,jj,ii)+c_1*ortho_1+c_2*ortho_2;
        end
        
        %update Locations (fast method)
        location2(:,:,ii)=cumsum(tangents2(:,:,ii)*length_link,2); 
        
        Xend2(ss,ii) = location2(1,end,ii);
        Yend2(ss,ii) = location2(2,end,ii);
    end
    comp_time=toc;
end
Xend2(Xend2==0)=NaN;
Yend2(Yend2==0)=NaN;

%signaling computation is finished
fprintf('\n>%u Configurations  with varaible iteration number completed, Computation 2 finished',P)
%% plots
   figure
for pp=1:P
    subplot(2,P,pp)
    Hx = hist(Xend(pp,:),bins(pp,:));
    bar(bins(pp,:),Hx);
    title(sprintf('[x] %i iterations and %i segments',defN,K(pp)));xlabel('position [nm]');ylabel('count')
    subplot(2,P,pp+P)
    Hy = hist(Yend(pp,:),bins(pp,:));hold on;
    bar(bins(pp,:),Hy);
    title(sprintf('[y] %i iterations and %i segments',defN,K(pp)));xlabel('position [nm]');ylabel('count') 

    %calculate the standard deviations of respective coordinates for each value of K
    devX(pp)=std(Xend(pp,:));
    devY(pp)=std(Yend(pp,:));
end

figure
for ss=1:S
    subplot(2,S,ss)
    Hx2 = hist(Xend2(ss,:),bins2);
    bar(bins2,Hx2);
    title(sprintf('[x] %i iterations and %i segments',N(ss),defK));xlabel('position [nm]');ylabel('count')
    subplot(2,S,ss+S)
    Hy2 = hist(Yend2(ss,:),bins2);hold on;
    bar(bins2,Hy2);
    title(sprintf('[y] %i iterations and %i segments',N(ss),defK));xlabel('position [nm]');ylabel('count') 

    devX2(ss) = std(Xend2(ss,1:N(ss)));
    devY2(ss) = std(Yend2(ss,1:N(ss)));
end
figure
subplot(2,2,1)
bar(bins(end,:),Hx);title(sprintf('zoom into K=%i segments (x)',K(end)))
xlabel('x [nm]')
subplot(2,2,2)
bar(bins(end,:),Hy);title(sprintf('zoom into K=%i segments (y)',K(end)))
xlabel('x [nm]')
subplot(2,2,3)
bar(bins2,Hx2);title(sprintf('zoom into N=%i iterations (x)',N(end)));
xlabel('x [nm]')
subplot(2,2,4)
bar(bins2,Hy2);title(sprintf('zoom into N=%i iterations (y)',N(end)));
xlabel('x [nm]')

figure
subplot(2,2,1)
plot(K,devX,K,sigma);title('standard deviation in x as a function of segment length');grid on;
legend('simulated deviation','predicted deviation');xlabel('number of segments')
subplot(2,2,2)
plot(K,devY,K,sigma);title('standard deviation in y as a function of segment length');grid on;
legend('simulated deviation','predicted deviation');xlabel('number of segments')
subplot(2,2,3)
plot(N,devX2,N,sigma2*ones(S,1));title('standard deviation in x as a function of iteration number');grid on;
legend('simulated deviation','predicted deviation');xlabel('iteration number')
subplot(2,2,4)
plot(N,devY2,N,sigma2*ones(S,1));title('standard deviation in y as a function of iteration number');grid on;
legend('simulated deviation','predicted deviation');xlabel('iteration number')

%closing statement
fprintf('\n>>>[Task 4] Completed.\n')
