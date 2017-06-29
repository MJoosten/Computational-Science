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
enable_debug_plots=true;
N=10000; %Iterations of Polymer/chain (DNA) generation (default:100)
K=2000; % Number of segments of chain (base pairs) (default:2000)
length_link=0.311;%[nm] Length of each chain link(base pair)(default:0.311)
length_persist=50; %[nm] persistence length (default:50)
length_chain=K*length_link; %[nm] Total length of chain (DNA)
t_initial=[1;0]; %initial orientation of t vector (unit length);
                 %Dont change this vector and expect this to work

%Preallocation
comp_time=zeros(N,2); %this array will hold the computational time for
                      %for each Iteration (N iterations)
distances=zeros(N,1); %will hold the squared end-to-end distances
location=zeros(2,K,N); %will hold the location for each polymer link

tangents=ones(2,K,N);% holds the angles %TODO: remove
average_distances = zeros(N,1);
% Calculate models 

%opening statement (for console iterpretability)
fprintf('\n>>>[task 1] Starting Computation with %u iterations and %u segments',N,K)


%generate random bend angles - mu=0;var=length_link/length_persistence
rand_angles=sqrt(length_link/length_persist)*randn(K,N);
angles_cum=cumsum(rand_angles);
cos_test=cos(angles_cum)';
sin_test=sin(angles_cum)';

for ii=1:N %loop over N iterations(generate N independent runs)
    tic %start a time for each run %FIX   

    tangents(:,:,ii)=[cos_test(ii,:);sin_test(ii,:)];
    
    comp_time(ii,1)=toc;    

    location(:,:,ii)=cumsum(tangents(:,:,ii)*length_link,2); 

    distances(ii)=sum(((location(:,end,ii)-location(:,1,ii))).^2);   
    comp_time(ii,2)=toc;
    average_distances(ii) = mean(distances(1:ii));
end

% Testing Stuff ----> remove please
predict_distance=4*length_persist*length_chain-8*length_persist^2*(1-exp(-length_chain/(2*length_persist)));
error=abs(predict_distance-mean(distances));
%%
if enable_plots %if you want to plot the generated WLC
%     close all
    figure % plot figures
    subplot(1,2,1) % all WLC chains
    for nn=1:min([N 100])
        scatter(location(1,:,nn),location(2,:,nn));hold on;
    end
    title(sprintf('WLC visualization of all %i chains',min([N 100])));xlabel('X position [nm]');ylabel('Y position [nm]')
    
    subplot(1,2,2) % sigle WL chain
    scatter(location(1,:,1),location(2,:,1),[],linspace(1,K,K),'filled')  
    title('WLC visualisation plot for the first iteration')
    xlabel('X position [nm]')
    ylabel('Y position [nm]')
end
if enable_debug_plots
    figure %debug plots
    subplot(2,2,1)
    plot(distances)
   title('Debugging distances - for first iteration of N')
   xlabel('Iterations (1:N)')
   ylabel('Square end to end distance')
    
    subplot(2,2,2) %testing biases in angles
   histogram(rand_angles(:,1))
   title('Debugging angles - for first iteration of N')
   xlabel('angle')
   ylabel('frequency')
   
   subplot(2,2,3)
   plot(1:N, average_distances,'color','green')
   refline(0,predict_distance)
   legend('computed average distances','expected average distances');
   title('average distances per iteration');xlabel('iteration number');
   ylabel('average square end-to-end distance');
end

% Printing Results --------------------------------------------------------

%print squared-end-to-end distances (mean + standard deviation)
%TODO: check if standard deviation is correct (see document)
fprintf('\n> The mean end-to-end distances for %u iterations is: %f, with standard deviation: %f',N,mean(distances),std(distances)/sqrt(N)) %

%print computational times
%TODO: figure out proper display values for times (comp times)
fprintf('\n> Total computational time for %u iterations is: %f, the mean time per iteration is:  %f',N,sum(comp_time(:,2)),mean(comp_time(:,2)));

%closing statement (for console iterpretability)
fprintf('\n>>> %u iterations completed, Computation finished\n',N)