%% Computational Science Final Project: Worm-Like Chain
% Task 4
% Authors: Maarten Joosten & Nemo Andrea
% IDs: xxxxxxx & 4473035
% Date of Creation: 20-06-2017
% github: https://github.com/MJoosten/Computational-Science

%% Task 3 - Prepping

clear all;
close all;
format compact;

%% Start

%parameters
enable_plots=true; %do you wish to plot the WLC? 
enable_debug_plots=false; %do you wish to view debugging plots? 
enable_error_plot=true;
plot_scatter=true; %do you wish to havea  scatter plot over plot3 for 
                   %the plot with multiple WLC in one figure?
P=8; %number of configurations (configuration = amount segments of chain) 
P_range=[50,10000]; %range of segment numbers (default: 100,5000)
N=1000; %Iterations of Polymer/chain (DNA) generation (default:100)
K=round(linspace(P_range(1),P_range(2),P)); % Number of segments of chain
                                     %(base pairs) (default:2000)
length_link=0.311;%[nm] Length of each chain link(base pair)(default:0.311)
length_persist=50; %[nm] persistence length (default:50)
length_chain=K*length_link; %[nm] Total length of chain (DNA)
t_initial=[0;0;1]; %initial orientation of t vector (unit length);
                   %(default: 0,0,1 (z axis))

                   
%Preallocation - Outside Loop
comp_time=zeros(N,1); %this array will hold the computational time for
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
    tangents=ones(3,K_local,N);% holds the angles %TODO: do this more efficiently

    % generate random bend angles
    % Gaussian Distribution with mu=0;var=length_link/length_persistence
    rand_angles1=sqrt(length_link/length_persist)*randn(K_local,N);
    rand_angles2=sqrt(length_link/length_persist)*randn(K_local,N);  
    %Cumulative Angles (these are the random angles values we will rotate at)
    cum_angles1=cumsum(rand_angles1,1);
    cum_angles2=cumsum(rand_angles2,1);
    
    %compute cosines and sines beforehand
    cos_test_1=cos(cum_angles1);
    sin_test_1=sin(cum_angles1);
    cos_test_2=cos(cum_angles2);
    sin_test_2=sin(cum_angles2);   
   

    % Computation ------------------------------------------------------------- 

    fprintf('\nComputing WLC 3D Distance for K=%u links, for N=%u iterations',K_local,N)
    for ii=1:N %loop over N iterations(generate N independent runs)
        tic %start a clock for each run        
             
        %compute tangents
        factor_cont=(sqrt(1-(sin_test_1(:,ii).*sin_test_2(:,ii)).^2));
        tangents(1,:,ii)=sin_test_1(:,ii).*cos_test_2(:,ii)./factor_cont;
        tangents(2,:,ii)=cos_test_1(:,ii).*sin_test_2(:,ii)./factor_cont;
        tangents(3,:,ii)=cos_test_1(:,ii).*cos_test_2(:,ii)./factor_cont;
       
        %update Locations (fast method)
        location(:,:,ii)=cumsum(tangents(:,:,ii)*length_link,2); 
        
        %Compute the squared end-to-end distance (works for non-zero starting
        %points too. Alternative method would be norm(vector)^2.      
        distances(ii,pp)=sum((location(:,end,ii)-location(:,1,ii)).^2);
        comp_time(ii)=toc; %clock in computation time for this (single) WLC.    
    end

    %calculate theoretical distance
    theoretical_dist=2*length_persist*length_chain-2*length_persist^2*...
        (1-exp(-length_chain/(length_persist)));
    theoretical_approx=2*length_persist*length_chain-2*length_persist^2;
    
    %calculate a percentage error using |new-old|/old 
    %<<<Keep in mind the figures only display the LAST iteration over P (last
    %configuration!>>>
    difference_percent=100*(abs(mean(distances(:,end))-...
        theoretical_dist(end))/theoretical_dist(end));
end


%signaling computation is finished
fprintf('\n>%u Configurations each with %u iterations completed, Computation finished',P,N)

% Displaying Results / PostProcessing -------------------------------------

%% Plotting 

%<<<Keep in mind the figures only display the LAST iteration over P (last
%configuration!>>>
if enable_debug_plots
    figure
    subplot(2,2,1)    
    histogram(rand_angles1,[-0.4 -0.4:0.02:0.4 0.4])
    title('[Task3]Distribution of angles (angle1)')
    xlabel('Angle')
    ylabel('Frequency')
    yyaxis right
   ylabel('probability')
   range_hist=-0.4:0.02:0.4;
   plot(range_hist,normpdf(range_hist,0,sqrt(length_link/length_persist)))
   legend('Computed Angles','Theoretical Distribution')
   
    subplot(2,2,2)
    histogram(rand_angles2,[-0.4 -0.4:0.02:0.4 0.4])
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
    %plotting a subset (multiple) WLC (set by plot_wlc_number)    
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
fprintf('\n> Total computational time for %u iterations is: %f, the mean time per iteration is:  %f',N,sum(comp_time),mean(comp_time));

%closing statement (for console iterpretability)
fprintf('\n>>>[task3] Completed.\n')