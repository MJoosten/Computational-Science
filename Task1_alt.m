%% Computational Science Final Project: Worm-Like Chain
% Task 1
% Authors: Maarten Joosten & Nemo Andrea
% IDs: xxxxxxx & 4473035
% Date of Creation: 20-06-2017
% github: https://github.com/MJoosten/Computational-Science

%% Task 1 -----------------------------------------------------------------

%prepping
clear all
close all
format compact

%% Start

enable_plots=true; %do you wish to plot the WLC? %debugging
N=100; %Iterations of Polymer/chain (DNA) generation (default:100)
K=2000; % Number of segments of chain (base pairs) (default:2000)
length_link=0.311;%[nm] Length of each chain link(base pair)(default:0.311)
length_persist=50; %[nm] persistence length (default:50)
length_chain=K*length_link; %[nm] Total length of chain (DNA)
t_initial=[1;0]; %initial orientation of t vector (unit length);

%Preallocation
comp_time=zeros(N,1); %this array will hold the computational time for
                      %for each Iteration (N iterations)
distances=zeros(N,1); %will hold the squared end-to-end distances
location=zeros(2,K,N); %will hold the location for each polymer link

tangents=ones(2,K,N);% holds the angles %TODO: remove
tangents(1,:,:)=tangents(1,:,:)*t_initial(1);
tangents(2,:,:)=tangents(2,:,:)*t_initial(2);% 

% Calculate models 

%opening statement (for console iterpretability)
fprintf('\n>>>[task 1] Starting Computation with %u iterations and %u segments',N,K)


%generate random bend angles - mu=0;var=length_link/length_persistence
rand_angles=sqrt(length_link/length_persist)*randn(N,K);

% % % angle_matrix=[cos(rand_angles),-sin(rand_angles);sin(rand_angles),cos(rand_angles)]
for ii=1:N %loop over N iterations(generate N independent runs)
    tic %start a time for each run %FIX
    for jj=1:K-1 %compute K segments
        tangents(:,jj+1,ii)=[cos(rand_angles(ii,jj)),-sin(rand_angles(ii,jj));sin(rand_angles(ii,jj)),cos(rand_angles(ii,jj))]*tangents(:,jj,ii);
        location(:,jj+1,ii)=location(:,jj,ii)+tangents(:,jj+1,ii)*length_link;
        
    end
    distances(ii)=sum((location(:,end,ii)-location(:,1,ii)).^2);   
    comp_time(ii)=toc;
    
end

if enable_plots %if you want to plot the generated WLC
    close all
    figure
    scatter(location(1,:,1),location(2,:,1),[],linspace(1,K,K),'filled')  
    title('WLC visualisation plot for the first iteration')
    xlabel('X position [nm]')
    ylabel('Y position [nm]')
    
    figure %testing biases between runs
    plot(distances)
   title('Debugging distances - for first iteration of N')
   xlabel('Iterations (1:N)')
   ylabel('Square end to end distance')
    
    figure %testing biases in angles
   histogram(rand_angles(1,:))
   title('Debugging angles - for first iteration of N')
   xlabel('angle')
   ylabel('frequency')
end

% Printing Results --------------------------------------------------------

%print squared-end-to-end distances (mean + standard deviation)
%TODO: check if standard deviation is correct (see document)
fprintf('\n> The mean end-to-end distances for %u iterations is: %f, with standard deviation: %f',N,mean(distances),std(distances/sqrt(N))) %

%print computational times
%TODO: figure out proper display values for times (comp times)
fprintf('\n> Total computational time for %u iterations is: %f, the mean time per iteration is:  %f',N,sum(comp_time),mean(comp_time));

%closing statement (for console iterpretability)
fprintf('\n>>> %u iterations completed, Computation finished\n',N)



% Testing Stuff ----> remove please
predict_distance=4*length_persist*length_chain-8*length_persist^2*(1-exp(-length_chain/(2*length_persist)))
error=abs(predict_distance-mean(distances))

