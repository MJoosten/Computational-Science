%% Computational Science Final Project: Worm-Like Chain
% Task 1
% Authors: Maarten Joosten & Nemo Andrea
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
Q=30; %how many different values for N do you wish to try?

N_range=[1,4]; %what range of values of N do you wish to try?
N=round(logspace(N_range(1),N_range(2),Q));
K=2000; % Number of segments of chain (base pairs) (default:2000)
length_link=0.311;%[nm] Length of each chain link(base pair)(default:0.311)
length_persist=50; %[nm] persistence length (default:50)
length_chain=K*length_link; %[nm] Total length of chain (DNA)
t_initial=[1;0]; %initial orientation of t vector (unit length);            
Q_avg=zeros(Q,1); %holds the average distance for each iteration of Q

fprintf('\n>>>[task 1] Starting Task 1')
for qq=1:Q
    
    N_local=N(qq);
    
    %Preallocation
    comp_time=0; %this will hold the computational time 
    distances=zeros(N_local,1); %will hold the squared end-to-end distances
    location=zeros(N_local,2,K); %will hold the location for each polymer link 
    tangents=ones(N_local,2,K);% holds the angles 
    average_distances = zeros(N_local,1);
   

    %opening statement (for console iterpretability)
    fprintf('\n> Starting Computation with %u iterations and %u segments',N_local,K)


    %generate random bend angles - mu=0;var=length_link/length_persistence
    rand_angles=sqrt(length_link/length_persist)*randn(K,N_local);

    %use cumulative approach to angles (cumulative rotation around z axis)
    angles_cum=cumsum(rand_angles);
    cos_test=cos(angles_cum)';
    sin_test=sin(angles_cum)';

    %generate rotation for starting x coordinate
    rotation_x=zeros(N_local,2,K);
    rotation_x(:,1,:)=cos_test;
    rotation_x(:,2,:)=sin_test;

    %generate rotation for starting y coordinate
    rotation_y=zeros(N_local,2,K);
    rotation_y(:,1,:)=-sin_test;
    rotation_y(:,2,:)=cos_test;

    %Doing the actual computation
    
    tic 
    fprintf('\n> Starting Vectorised Computation')

    tangents=rotation_x*t_initial(1)+rotation_y*t_initial(2);

    location=cumsum(tangents*length_link,3); 

    distances=sum(((location(:,:,1)-location(:,:,end))).^2,2);   

    comp_time=toc;
    fprintf('\n> Finished Vectorised Computation: %f seconds',comp_time)
    
    Q_avg(qq)=mean(distances);
end
% Debugging Calculations

%compute cumulative means
running_average_distances = cumsum(distances')./(1:length(distances));
predict_distance=4*length_persist*length_chain-8*length_persist^2*(1-exp(-length_chain/(2*length_persist)));
error=abs(predict_distance-mean(distances));
%% Plotting Section

if enable_plots %if you want to plot the generated WLC
     close all
    figure % plot figures
    subplot(1,2,1) % all WLC chains
    for nn=1:min([N(end) 100])
        scatter(location(nn,1,:),location(nn,2,:));hold on;
    end
    title(sprintf('[TASK 1]WLC visualization of all %i chains',min([N(end) 100])));xlabel('X position [nm]');ylabel('Y position [nm]')
    
    subplot(1,2,2) % sigle WL chain
    scatter(location(1,1,:),location(1,2,:),[],linspace(1,K,K),'filled')  
    title('[TASK 1]WLC visualisation plot for the first iteration')
    xlabel('X position [nm]')
    ylabel('Y position [nm]')
end
if enable_debug_plots
    figure %debug plots
    subplot(2,2,1)
    plot(distances)
   title('[TASK 1]Debugging distances - for first iteration of N')
   xlabel('Iterations (1:N)')
   ylabel('Square end to end distance')
    
    subplot(2,2,2) %testing biases in angles
   histogram(rand_angles(:,1),[-0.4 -0.4:0.02:0.4 0.4])
   title('[TASK 1]Debugging angles - for first iteration of N')
   xlabel('angle')
   yyaxis left
   ylabel('frequency')
   yyaxis right
   ylabel('probability')
   range_hist=-0.4:0.02:0.4;
   plot(range_hist,normpdf(range_hist,0,sqrt(length_link/length_persist)))
   legend('Computed Angles','Theoretical Distribution')
   
   subplot(2,2,3)
   plot(1:N(end), running_average_distances,'color','green')
   refline(0,predict_distance)
   legend('computed average distances','expected average distances');
   title('[TASK 1] Debugging average distances per iteration');xlabel('iteration number');
   ylabel('average square end-to-end distance');
   
   subplot(2,2,4)
   loglog(N,abs(Q_avg-predict_distance),N,1./sqrt(N))  
   grid on
   title('[TASK 1] Viewing 1/sqrt(N) dependence monte carlo')
   xlabel('Number of WLC distances averaged')
   ylabel('Absolute error from expected value')
   legend('Computed data', 'f(N)=1/sqrt(N))')
end

% Printing Results --------------------------------------------------------

%print squared-end-to-end distances (mean + standard deviation)
%TODO: check if standard deviation is correct (see document)
fprintf('\n> The mean end-to-end distances for %u iterations is: %f, with standard deviation: %f',N(end),mean(distances),std(distances)/sqrt(N(end))) %

%print computational times
%TODO: figure out proper display values for times (comp times)
fprintf('\n> Total computational time for %u iterations is: %f, the mean time per iteration is:  %f',N(end),comp_time,comp_time/N(end));

%closing statement (for console iterpretability)
fprintf('\n>>> %u different N values with max %u iterations completed, Computation finished. Task 1 Done. \n',Q,N(end))

