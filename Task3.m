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