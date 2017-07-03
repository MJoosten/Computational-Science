%% Computational Science Final Project: Worm-Like Chain
% Task 1
% Authors: Maarten Joosten & Nemo Andrea
% IDs: 4662661 & 4473035
% Date of Creation: 20-06-2017
% github: https://github.com/MJoosten/Computational-Science

%prepping
clear all
close all
format compact
%% Task 1 -----------------------------------------------------------------
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

%compute cumulative means
running_average_distances = cumsum(distances')./(1:length(distances));
predict_distance=4*length_persist*length_chain-8*length_persist^2*(1-exp(-length_chain/(2*length_persist)));
error=abs(predict_distance-mean(distances));
%% Plotting Section
if enable_plots %if you want to plot the generated WLC
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
%%
sprintf('task 1 finished, press enter to start next task')
pause()
clc
%% Task 2 -----------------------------------------------------------------
%% Start
enable_plots=true; %do you wish to plot the WLC? %debugging
enable_speed_plots=false; %do you wish to evaluate running times?
P=20; %number of different values for K to try (affects total length)
P_range=[2000,4000]; %(default:1500,20000)
N=250; %Iterations of Polymer/chain (DNA) generation (default:100)
K=round(linspace(P_range(1),P_range(2),P)); % Number of segments of chain (base pairs)
length_link=0.311;%[nm] Length of each chain link(base pair)(default:0.311)
length_persist=30; %[nm] persistence length (default:50)
length_chain=K*length_link; %[nm] Total length of chain (DNA)
t_initial=[1;0]; %initial orientation of t vector (unit length);
time_comp=zeros(P,1);
%Some Preallocation
distances=zeros(N,P);
distances_to_4=zeros(N,P);

% Calculate models -------------------------------------------------------- 
%opening statement (for console iterpretability)
fprintf(['\n>>>[task 2] Starting Computation with %u iterations,'...
            '%u Configurations (P) and %u to %u segments'],N,P,K(1),K(end))

for pp=1:P
    tic;
    K_select=K(pp); %select the value for K (#links)    
    
    %Main Preallocation
    location=zeros(N,2,K_select); %will hold the location for each polymer link
    tangents=ones(N,2,K_select);% holds the angles 
    
    %generate random bend angles - mu=0;var=length_link/length_persistence    
    rand_angles=sqrt(length_link/length_persist)*randn(K_select,N);
    angles_cum=cumsum(rand_angles);
    cos_test=cos(angles_cum)';
    sin_test=sin(angles_cum)';
    
    %generate rotation for starting x coordinate
    rotation_x=zeros(N,2,K_select);
    rotation_x(:,1,:)=cos_test;
    rotation_x(:,2,:)=sin_test;

    %generate rotation for starting y coordinate
    rotation_y=zeros(N,2,K_select);
    rotation_y(:,1,:)=-sin_test;
    rotation_y(:,2,:)=cos_test;

    fprintf('\nComputing WLC Distance for K=%u links, for N=%u iterations',K_select,N)
     
    tangents=rotation_x*t_initial(1)+rotation_y*t_initial(2);

    location=cumsum(tangents*length_link,3); 

    distances(:,pp)=sum(((location(:,:,1)-location(:,:,end))).^2,2);  
    
    time_comp(pp)=toc;
end
%closing statement (for console iterpretability)
fprintf('\n> %u configurations completed, Computation finished\n',P)

%compute theoretical values
theoretical_full=4*length_persist*length_chain-8*length_persist^2*...
    (1-exp(-length_chain/(length_persist*2)));    
theoretical_approx=4*length_persist*length_chain-8*length_persist^2;
for ii=1:P    
    error_chain(ii)=sqrt((mean(distances(:,ii).^2)-mean(distances(:,ii))^2)/N);
end

%% plotting figures
if enable_plots
    figure
    subplot(1,2,1)
    errorbar(length_chain,mean(distances,1),error_chain);
    hold on
    plot(length_chain,theoretical_full)
    hold on
    plot(length_chain,theoretical_approx)
    hold off
    title('[Task 2]Length of chain versus the squared end-to-end distance')
    xlabel('Length of chain [nm]')
    ylabel('squared end to end distance [nm^2]')
    legend('Monte Carlo','Theoretical Values','Theoretical Values (using approximation)')

    subplot(1,2,2)
    plot(length_chain,exp(-length_chain./(length_persist*2)))
    title('f(L) = exp(-L/2*ksi)');xlabel('Length of chain [nm]');ylabel('f(L)');
end

%to evaluate the computational time per link. We expect a linear
%relationship
if enable_speed_plots
    figure
    plot(K,time_comp)
    title('[Task 2]computational time per link count')
    xlabel('number of links')
    ylabel(sprintf('time to compute WLC (2D) for %u iterations at given K [seconds]',N))
end
%% Fitting Parameter

%end to end distance
distances_end=mean(distances,1);
%estimate the possible range for persistence length
unknown=linspace(1,500,200);

for ii=1:200
% %     Chi_Deriv(ii)=sum(2*(16*unknown(ii)-4*length_chain).*...
% %         (mean(distances,1)-4*unknown(ii).*length_chain+8*unknown(ii).^2));
Chi_Deriv(ii)=Chi_squared_deriv(unknown(ii),length_chain,distances_end,error_chain);
end
figure

plot(unknown,Chi_Deriv,'red')
hold on
xlabel('Persistence Length [nm]')
ylabel('\chi^2')
refline(0,0)

title(sprintf('[Task 2]Chi^2 as as function of persistence length; true value: %i',length_persist))


estimate_ksi=unknown(1); %start at max of user provided range
error=10;
tic
counter=1;
while (error>10^-3&&toc<5)
    estimate_ksi=estimate_ksi-Chi_squared_deriv(estimate_ksi,length_chain,distances_end,error_chain)/(Chi_second_deriv(estimate_ksi,length_chain,distances_end,error_chain));
    error=abs(Chi_squared_deriv(estimate_ksi,length_chain,distances_end,error_chain)); %deviation from 0 (abs(y))%update a new error 
    estimate_array(counter,1)=estimate_ksi;
    estimate_array(counter,2)=Chi_squared_deriv(estimate_ksi,length_chain,distances_end,error_chain);
    counter=counter+1; %keep track of how many iterations things took
end
plot(estimate_array(:,1),estimate_array(:,2),'c*')
fprintf('> Root Solving using Newton-Raphson Method completed in %i iterations',counter)
fprintf('> Estimate for persistence length is: %f nm. True value was: %i nm.',estimate_ksi,length_persist)
line([estimate_ksi,estimate_ksi],ylim)
legend('Chi Squared','Zero Line','Root Solving Points','Estimate')

fprintf('\n>>>[Task 2] completed.\n',P)
%%
sprintf('task 2 finished, press enter to start next task')
pause()
clc
%% Task 3 - Prepping
%% Start
%----------------------TODO: VECTORISE N, BABY ----------------------------
%parameters
enable_plots=true; %do you wish to plot the WLC? 
enable_debug_plots=true; %do you wish to view debugging plots? 
enable_error_plot=true;
plot_scatter=true; %do you wish to havea  scatter plot over plot3 for 
                   %the plot with multiple WLC in one figure?
                   
                  
P=50; %number of configurations (configuration = amount segments of chain) 
P_range=[10,1000]; %range of segment numbers (default: 100,5000)
N=500; %Iterations of Polymer/chain (DNA) generation (default:100)
K=round(linspace(P_range(1),P_range(2),P)); % Number of segments of chain
                                     %(base pairs) (default:2000)
length_link=0.311;%[nm] Length of each chain link(base pair)(default:0.311)
length_persist=35; %[nm] persistence length (default:50)
length_chain=K*length_link; %[nm] Total length of chain (DNA)
t_initial=[0;0;1]; %initial orientation of t vector (unit length);
                   %(default: 0,0,1 (z axis))

                   
%Preallocation - Outside Loop
comp_time=0; %this array will hold the computational time for
                      %for each Iteration (N iterations)
distances=zeros(N,P); %will hold the squared end-to-end distances

%opening statement (for console iterpretability)
fprintf(['\n>>>[task 3] Starting Computation WLC 3D with %u'...
        ' configurations each with %u iterations and number of segments '...
        'between %u and %u'],P,N,min(K),max(K))
for pp=1:P
    K_local=K(pp);
    
    %Preallocation - Inside Loop
    location=zeros(3,K_local,N); %will hold the location for each polymer link (3D)
    tangents=ones(3,K_local,N);% holds the angles 
    norm_factor=zeros(N,1);
    ortho_1=zeros(N,1);
    ortho_2=zeros(N,1);
    alpha_t=zeros(N,1);
    beta_t=zeros(N,1);
    c_t=zeros(N,1);
    c_1=zeros(N,1);
    c_2=zeros(N,1);
    
    %TODO: do this more efficiently
    tangents(1,:,:)=tangents(1,:,:)*t_initial(1); %setting initial tangent
    tangents(2,:,:)=tangents(2,:,:)*t_initial(2); %setting initial tangent
    tangents(3,:,:)=tangents(3,:,:)*t_initial(3); %setting initial tangent

    
    % generate random bend angles
    % Gaussian Distribution with mu=0;var=length_link/length_persistence
     rand_angles=sqrt(length_link/length_persist)*randn(2,K_local,N);
     cos_1=reshape(cos(rand_angles(1,:,:)),[K_local,N]);
     sin_1=reshape(sin(rand_angles(1,:,:)),[K_local,N]);
     cos_2=reshape(cos(rand_angles(2,:,:)),[K_local,N]);
     sin_2=reshape(sin(rand_angles(2,:,:)),[K_local,N]);
     
    % Computation ------------------------------------------------------------- 

    fprintf('\nComputing WLC 3D Distance for K=%u links, for N=%u iterations',K_local,N)
    tic %start a clock for each run 
           
    
    for jj=1:K_local-1 %compute K segments %FIX                      
        %find alpha and beta of PREVIOUS iteration
        alpha_t=reshape(acos(tangents(3,jj,:)),[1,N]); %arccos(t_z)       
        beta_t=reshape(atan2(tangents(2,jj,:),tangents(1,jj,:)),[1,N]);%arctan(t_y/t_x)            
       
        ortho_1=[cos(alpha_t).*cos(beta_t);cos(alpha_t).*sin(beta_t);-sin(alpha_t)];
        ortho_2=[-sin(beta_t);cos(beta_t);zeros(1,N)];

        %calculate coefficients       
        norm_factor=sqrt(1-(sin_1(jj,:).*sin_2(jj,:)).^2);
        c_t=(cos_1(jj,:).*cos_2(jj,:))./norm_factor;
        c_1=(sin_1(jj,:).*cos_2(jj,:))./norm_factor;
        c_2=(cos_1(jj,:).*sin_2(jj,:))./norm_factor;        

        %calculate the new tangent vector (3D)
        tangents(1,jj+1,:)=c_t.*reshape(tangents(1,jj,:),[1,N])+c_1.*ortho_1(1,:)+c_2.*ortho_2(1,:);
        tangents(2,jj+1,:)=c_t.*reshape(tangents(2,jj,:),[1,N])+c_1.*ortho_1(2,:)+c_2.*ortho_2(2,:);
        tangents(3,jj+1,:)=c_t.*reshape(tangents(3,jj,:),[1,N])+c_1.*ortho_1(3,:)+c_2.*ortho_2(3,:);
    end

    %update Locations (fast method)
    location=cumsum(tangents*length_link,2); 

    %Compute the squared end-to-end distance (works for non-zero starting
    %points too. Alternative method would be norm(vector)^2.      
    distances(:,pp)=sum((location(:,end,:)-location(:,1,:)).^2);
           
    
    comp_time=toc; %clock in computation time for this  WLC set. 
    
    %calculate theoretical distance
    theoretical_dist=2*length_persist*length_chain-2*length_persist^2*...
        (1-exp(-length_chain/length_persist));
    theoretical_approx=2*length_persist*length_chain-2*length_persist^2;
    
    %calculate a percentage error using |new-old|/old 
    %<<<Keep in mind the figures only display the LAST iteration over P (last
    %configuration!>>>
    difference_percent=100*(abs(mean(distances(:,end))-...
        theoretical_dist(end))/theoretical_dist(end));
end
%signaling computation is finished
fprintf('\n>%u Configurations each with %u iterations completed, Computation finished',P,N)

%% Plotting 

%<<<Keep in mind the figures only display the LAST iteration over P (last
%configuration!>>>
if enable_debug_plots
    figure
    subplot(2,2,1)    
    histogram(rand_angles(1,:,:),[-0.4 -0.4:0.02:0.4 0.4])
    title('[Task3]Distribution of angles (angle1)')
    xlabel('Angle')
    ylabel('Frequency')
    yyaxis right
   ylabel('probability')
   range_hist=-0.4:0.02:0.4;
   plot(range_hist,normpdf(range_hist,0,sqrt(length_link/length_persist)))
   legend('Computed Angles','Theoretical Distribution')
   
    subplot(2,2,2)
    histogram(rand_angles(2,:,:),[-0.4 -0.4:0.02:0.4 0.4])
    title('[Task3]Distribution of angles (angle2)')
    xlabel('Angle')
    yyaxis right
   ylabel('probability')
   range_hist=-0.4:0.02:0.4;
   plot(range_hist,normpdf(range_hist,0,sqrt(length_link/length_persist)))
   legend('Computed Angles','Theoretical Distribution')
    
    ylabel('Frequency')
    subplot(2,2,3)
    plot(distances)
    title('[Task3]Distances for different WLC (no bias expected)')
    xlabel('Iteration number (N)')
    ylabel('Squared end-to-end distance [nm^2]')
    
    %TODO: add curvature??? or More debuggin options
end

%Plotting Results

%<<<Keep in mind the figures only display the LAST iteration over P (last
%configuration!>>>
if enable_plots
    figure
    %plotting a subset (multiple) WLC  
    subplot(1,2,1)
    if plot_scatter
        for ii=1:min([N 100])
            scatter3(location(1,:,ii),location(2,:,ii),location(3,:,ii),...
                [],linspace(1,K(end),K(end)),'filled')
            hold on
        end
        title(sprintf('[Task 3]WLC plots for the first %u iterations',min([N 100])))
        xlabel('X position [nm]')
        ylabel('Y position [nm]')
        zlabel('Z position [nm]')
    else
        for ii=1:min([N 100])
            plot3(location(1,:,ii),location(2,:,ii),location(3,:,ii))
            hold on
        end
        title(sprintf('[Task 3]WLC plots for the first %u iterations',min([N 100])))
        grid on
        xlabel('X position [nm]')
        ylabel('Y position [nm]')
        zlabel('Z position [nm]')
    end
    %plotting a single WLC 
    subplot(1,2,2)
    scatter3(location(1,:,1),location(2,:,1),location(3,:,1),...
        [],linspace(1,K(end),K(end)),'filled')
    title('[Task 3]WLC plot for the first iteration (single WLC)')
    xlabel('X position [nm]')
    ylabel('Y position [nm]')
    zlabel('Z position [nm]')
end

%Plotting Error plot (3B)

if enable_error_plot
    %%make non-for loop? more elegant, no need for performance
    %%TODO: can you check this? it seems okay now, but id like to have someone
  
    clear error_chain
    for ii=1:P    
        error_chain(ii)=sqrt((mean(distances(:,ii).^2)-mean(distances(:,ii))^2)/N);
        %no need for preallocation here
    end   
   
    figure
    errorbar(length_chain,mean(distances,1),error_chain);
    hold on
    plot(length_chain,theoretical_dist)
    hold on
    plot(length_chain,theoretical_approx)
    hold off
    title(sprintf('[Task 3]Length of chain versus the squared end-to-end distance [3D]; N=%u',N))
    xlabel('Length of chain [nm]')
    ylabel('squared end to end distance [nm^2]')
    legend('Monte Carlo','Theoretical Values','Theoretical Values (using approximation)')
end

%Printing Some Messages (Task 3A) -----------------------------------------
fprintf(['\nNote: The following statements only apply to the last '...
         'iteration over P, that is K=%u'],K(end))
%print squared-end-to-end distances (mean + standard deviation)
%TODO: check if standard deviation is correct (see document)
%TODO: check how we wish to display values such as percentage (how many
%decimals etc)
%<<<Keep in mind the figures only display the LAST iteration over P (last
%configuration!>>>
fprintf('\n> The mean squared end-to-end distances for %u iterations is: %f, with standard deviation: %f',N,mean(distances(:,end)),std(distances(:,end))/sqrt(N)) %
fprintf('\n> The theoretical squared end-to-end distance is: %f; the difference in percentage is: %u',theoretical_dist(end),round(difference_percent));

%print computational times
%TODO: figure out proper display values for times (comp times)
%<<<Keep in mind the figures only display the LAST iteration over P (last
%configuration!>>>
fprintf('\n> Total computational time for %u iterations is: %f, the mean time per iteration is:  %f',N,comp_time,comp_time/N);

%% Fitting Parameter

%end to end distance
distances_end=mean(distances,1);
%estimate the possible range for persistence length
unknown=linspace(1,500,200);

for ii=1:200
% %     Chi_Deriv(ii)=sum(2*(16*unknown(ii)-4*length_chain).*...
% %         (mean(distances,1)-4*unknown(ii).*length_chain+8*unknown(ii).^2));
Chi_Deriv(ii)=Chi_squared_deriv(unknown(ii),length_chain,distances_end,error_chain);
end
figure

plot(unknown,Chi_Deriv,'red')
hold on
xlabel('Persistence Length [nm]')
ylabel('\chi^2')
refline(0,0)

title(sprintf('[Task 3] Chi^2 as as function of persistence length; true value: %i',length_persist))


estimate_ksi=unknown(1); %start at max of user provided range
error=10;
tic
counter=1;
while (error>10^-3&&toc<5)
    estimate_ksi=estimate_ksi-Chi_squared_deriv(estimate_ksi,length_chain,distances_end,error_chain)/(Chi_second_deriv(estimate_ksi,length_chain,distances_end,error_chain));
    error=abs(Chi_squared_deriv(estimate_ksi,length_chain,distances_end,error_chain)); %deviation from 0 (abs(y))%update a new error 
    estimate_array(counter,1)=estimate_ksi;
    estimate_array(counter,2)=Chi_squared_deriv(estimate_ksi,length_chain,distances_end,error_chain);
    counter=counter+1; %keep track of how many iterations things took
end
plot(estimate_array(:,1),estimate_array(:,2),'c*')
fprintf('\n> Root Solving using Newton-Raphson Method completed in %i iterations',counter)
fprintf('\n> Estimate for persistence length is: %f nm. True value was: %i nm.',estimate_ksi,length_persist)
line([estimate_ksi,estimate_ksi],ylim)
legend('Chi Squared','Zero Line','Root Solving Points','Estimate')

%closing statement (for console iterpretability)
fprintf('\n>>>[task3] Completed.\n')
%%
sprintf('task 3 finished, press enter to start next task')
pause()
clc
%% Task 4
% Start

%parameters
enable_plots=true; %do you wish to plot the WLC? %debugging
enable_debug_plots=true;
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
if enable_debug_plots==true
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
end

if enable_plots==true
figure
subplot(2,2,1)
bar(bins(end,:),Hx);title(sprintf('[task 4] zoom into K=%i segments (x)',K(end)))
xlabel('x [nm]')
subplot(2,2,2)
bar(bins(end,:),Hy);title(sprintf('[task 4] zoom into K=%i segments (y)',K(end)))
xlabel('x [nm]')
subplot(2,2,3)
bar(bins2,Hx2);title(sprintf('[task 4] zoom into N=%i iterations (x)',N(end)));
xlabel('x [nm]')
subplot(2,2,4)
bar(bins2,Hy2);title(sprintf('[task 4] zoom into N=%i iterations (y)',N(end)));
xlabel('x [nm]')

figure
subplot(2,2,1)
plot(K,devX,K,sigma);title('[task 4]standard deviation in x as a function of segment length');grid on;
legend('simulated deviation','predicted deviation');xlabel('number of segments')
subplot(2,2,2)
plot(K,devY,K,sigma);title('[task 4]standard deviation in y as a function of segment length');grid on;
legend('simulated deviation','predicted deviation');xlabel('number of segments')
subplot(2,2,3)
plot(N,devX2,N,sigma2*ones(S,1));title('[task 4]standard deviation in x as a function of iteration number');grid on;
legend('simulated deviation','predicted deviation');xlabel('iteration number')
subplot(2,2,4)
plot(N,devY2,N,sigma2*ones(S,1));title('[task 4]standard deviation in y as a function of iteration number');grid on;
legend('simulated deviation','predicted deviation');xlabel('iteration number')
end
%closing statement
fprintf('\n>>>[Task 4] Completed.\n')

%%
function y_out=Chi_squared_deriv(input,length_chain,distances,error_chain)
%find out how many iterations of P we had (output=P)
length_P=length(distances);
    for ii=1:length_P
        factor=1/(error_chain(ii)^2);
        term1=16*input^3*factor;
        term2=-24*length_chain(ii)*input^2*factor;
        term3=(8*length_chain(ii)^2+8*distances(ii))*input*factor;
        term4=-(4*length_chain(ii)*distances(ii))*factor;
    end
y_out=(term1+term2+term3+term4)/length_P;

end

function y_out=Chi_second_deriv(input,length_chain,distances,error_chain)
 length_P=length(distances);
    for ii=1:length_P
        factor=1/(error_chain(ii)^2);
        term1=3*16*input^2*factor;
        term2=-2*24*length_chain(ii)*input*factor;
        term3=(8*length_chain(ii)^2+8*distances(ii))*factor;
    % % %     term4=-(8*length_chain(ii)*distances(ii))*factor;
    end
y_out=term1+term2+term3;
end