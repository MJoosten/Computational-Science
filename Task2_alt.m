%% Computational Science Final Project: Worm-Like Chain
% Task 2
% Authors: Maarten Joosten & Nemo Andrea
% Date of Creation: 20-06-2017
% github: https://github.com/MJoosten/Computational-Science

%% Task 2 -----------------------------------------------------------------

%prepping
clear all
close all
format compact



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


%%make non-for loop? more elegant, no need for performance
%%TODO: can you check this? it seems okay now, but id like to have someone
%%check this
for ii=1:P    
    error_chain(ii)=sqrt((mean(distances(:,ii).^2)-mean(distances(:,ii))^2)/N);
end

% plotting figures
%close all
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
function y_out=Chi_squared_deriv(input,length_chain,distances,error_chain)
% % % y_out=sum(2*(16*input-4*length_chain).*(mean(distances)-4*input.*length_chain+8*input.^2));
% % % y_out=(256*input^3-...
% % %     192*sum(length_chain)*input^2+...
% % %     (32*sum(length_chain.^2)+32*sum(distances))*input-...
% % %     2*4*sum(length_chain.*distances))*...
% % %     sum(1./(error_chain.^2));
%find out how many iterations of P we had (output=P)
length_P=length(distances);
    for ii=1:length_P
        factor=1/(error_chain(ii)^2);
        term1=256*input^3*factor;
        term2=-192*length_chain(ii)*input^2*factor;
        term3=(32*length_chain(ii)^2+32*distances(ii))*input*factor;
        term4=-(8*length_chain(ii)*distances(ii))*factor;
    end
y_out=(term1+term2+term3+term4)/length_P;

end

function y_out=Chi_second_deriv(input,length_chain,distances,error_chain)
 length_P=length(distances);
    for ii=1:length_P
        factor=1/(error_chain(ii)^2);
        term1=3*256*input^2*factor;
        term2=-2*192*length_chain(ii)*input*factor;
        term3=(32*length_chain(ii)^2+32*distances(ii))*factor;
    % % %     term4=-(8*length_chain(ii)*distances(ii))*factor;
    end
y_out=term1+term2+term3;
end