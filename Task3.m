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

%----------------------TODO: VECTORISE, BABY ------------------------------
%parameters
enable_plots=true; %do you wish to plot the WLC? 
enable_debug_plots=true; %do you wish to view debugging plots? 
plot_scatter=false; %do you wish to havea  scatter plot over plot3 for 
                   %the plot with multiple WLC in one figure?
plot_wlc_number=10; %how many WLC do you wish to plot in one single figure?
                    %requires enable_plots to be true to have effect
N=100; %Iterations of Polymer/chain (DNA) generation (default:100)
K=2000; % Number of segments of chain (base pairs) (default:2000)
length_link=0.311;%[nm] Length of each chain link(base pair)(default:0.311)
length_persist=50; %[nm] persistence length (default:50)
length_chain=K*length_link; %[nm] Total length of chain (DNA)
t_initial=[0;0;1]; %initial orientation of t vector (unit length);
                   %(default: 0,0,1 (z axis))

%Preallocation
comp_time=zeros(N,1); %this array will hold the computational time for
                      %for each Iteration (N iterations)
distances=zeros(N,1); %will hold the squared end-to-end distances
location=zeros(3,K,N); %will hold the location for each polymer link (3D)
tangents=ones(3,K,N);% holds the angles %TODO: do this more efficiently
tangents(1,:,:)=tangents(1,:,:)*t_initial(1); %setting initial tangent
tangents(2,:,:)=tangents(2,:,:)*t_initial(2); %setting initial tangent
tangents(3,:,:)=tangents(3,:,:)*t_initial(3); %setting initial tangent

% generate random bend angles
% Gaussian Distribution with mu=0;var=length_link/length_persistence
rand_angles=sqrt(length_link/length_persist)*randn(2,K,N);

% Computation -------------------------------------------------------------

%opening statement (for console iterpretability)
fprintf('\n>>>[task 3] Starting Computation WLC 3D with %u iterations and %u segments',N,K)

for ii=1:N %loop over N iterations(generate N independent runs)
    tic %start a clock for each run 
    for jj=1:K-1 %compute K segments %FIX                      
        %find alpha and beta of PREVIOUS iteration
        alpha_t=acos(tangents(3,jj,ii)); %arccos(t_z)
        %beta_t=atan2(tangents(2,jj,ii)/tangents(1,jj,ii));
        %%arctan(t_y/t_x) %OLD
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
        
        %update the location of the WLC (3D)
        location(:,jj+1,ii)=location(:,jj,ii)+tangents(:,jj+1,ii)*length_link;        
    end
    %Compute the squared end-to-end distance (works for non-zero starting
    %points too. Alternative method would be norm(vector)^2.
    distances(ii)=sum((location(:,end,ii)-location(:,1,ii)).^2);   
    comp_time(ii)=toc; %clock in computation time for this (single) WLC.    
end

%calculate theoretical distance
theoretical_dist=2*length_persist*length_chain-2*length_persist^2*...
    (1-exp(-length_chain/length_persist));
%calculate a percentage error using |new-old|/old
difference_percent=100*(abs(mean(distances)-theoretical_dist)/...
    theoretical_dist);

%signaling computation is finished
fprintf('\n> %u iterations completed, Computation finished',N)

% Displaying Results / PostProcessing -------------------------------------

%Plotting Debugging plots
if enable_debug_plots
    figure
    subplot(2,2,1)    
    histogram(rand_angles(1,:,:))
    title('Distribution of angles (angle1)')
    xlabel('Angle')
    ylabel('Frequency')
    subplot(2,2,2)
    histogram(rand_angles(2,:,:))
    title('Distribution of angles (angle2)')
    xlabel('Angle')
    ylabel('Frequency')
    subplot(2,2,3)
    plot(distances)
    title('Distances for different WLC (no bias expected)')
    xlabel('Iteration number (N)')
    ylabel('Squared end-to-end distance [nm^2]')
    
    %TODO: add curvature??? or More debuggin options
end

%Plotting Results
if enable_plots
    figure
    %plotting a subset (multiple) WLC (set by plot_wlc_number)    
    subplot(1,2,1)
    if plot_scatter
        for ii=1:plot_wlc_number
            scatter3(location(1,:,ii),location(2,:,ii),location(3,:,ii),[],linspace(1,K,K),'filled')
            hold on
        end
        title(sprintf('WLC plots for the first %u iterations',plot_wlc_number))
        xlabel('X position [nm]')
        ylabel('Y position [nm]')
        zlabel('Z position [nm]')
    else
        for ii=1:plot_wlc_number
            plot3(location(1,:,ii),location(2,:,ii),location(3,:,ii))
            hold on
        end
        title(sprintf('WLC plots for the first %u iterations',plot_wlc_number))
        grid on
        xlabel('X position [nm]')
        ylabel('Y position [nm]')
        zlabel('Z position [nm]')
    end
    %plotting a single WLC 
    subplot(1,2,2)
    scatter3(location(1,:,1),location(2,:,1),location(3,:,1),[],linspace(1,K,K),'filled')
    title('WLC plot for the first iteration (single WLC)')
    xlabel('X position [nm]')
    ylabel('Y position [nm]')
    zlabel('Z position [nm]')
end

%print squared-end-to-end distances (mean + standard deviation)
%TODO: check if standard deviation is correct (see document)
%TODO: check how we wish to display values such as percentage (how many
%decimals etc)
fprintf('\n> The mean squared end-to-end distances for %u iterations is: %f, with standard deviation: %f',N,mean(distances),std(distances)/sqrt(N)) %
fprintf('\n> The theoretical squared end-to-end distance is: %f; the difference in percentage is: %u',theoretical_dist,round(difference_percent));

%print computational times
%TODO: figure out proper display values for times (comp times)
fprintf('\n> Total computational time for %u iterations is: %f, the mean time per iteration is:  %f',N,sum(comp_time),mean(comp_time));

%closing statement (for console iterpretability)
fprintf('\n>>>[task3] Completed.\n')
