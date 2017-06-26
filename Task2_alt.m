%% Computational Science Final Project: Worm-Like Chain
% Task 2
% Authors: Maarten Joosten & Nemo Andrea
% IDs: xxxxxxx & 4473035
% Date of Creation: 20-06-2017
% github: https://github.com/MJoosten/Computational-Science

%% Task 2 -----------------------------------------------------------------

%prepping
clear all
close all
format compact



%% Start

enable_plots=true; %do you wish to plot the WLC? %debugging
P=8; %number of different values for K to try (affects total length)
P_range=[1500,20000]; %(default:1500,20000)
N=100; %Iterations of Polymer/chain (DNA) generation (default:100)
K=round(linspace(P_range(1),P_range(2),P)); % Number of segments of chain (base pairs)
length_link=0.311;%[nm] Length of each chain link(base pair)(default:0.311)
length_persist=85; %[nm] persistence length (default:50)
length_chain=K*length_link; %[nm] Total length of chain (DNA)
t_initial=[1;0]; %initial orientation of t vector (unit length);
time=0;
%Some Preallocation
distances=zeros(N,P);
distances_to_4=zeros(N,P);

% Calculate models -------------------------------------------------------- 

%opening statement (for console iterpretability)
fprintf('\n>>>[task 2] Starting Computation with %u iterations, %u Configurations (P) and %u to %u segments',N,P,K(1),K(end))

for pp=1:P
    K_select=K(pp);
    
    %Main Preallocation
    location=zeros(2,K_select,N); %will hold the location for each polymer link
    tangents=ones(2,K_select,N);% holds the angles %TODO: remove
    tangents(1,:,:)=tangents(1,:,:)*t_initial(1);
    tangents(2,:,:)=tangents(2,:,:)*t_initial(2);% 
    rand_angles=sqrt(length_link/length_persist)*randn(K_select,N);
    cos_t=cos(rand_angles);
    sin_t=sin(rand_angles);    
    
    fprintf('\nComputing WLC Distance for K=%u links, for N=%u iterations',K_select,N)
    tic;
    for ii=1:N %loop over N iterations(generate N independent runs)       
        for jj=1:K_select-1 %compute K segments
            tangents(:,jj+1,ii)=[cos_t(jj,ii),-sin_t(jj,ii);sin_t(jj,ii),cos_t(jj,ii)]*tangents(:,jj,ii);           
        end
        
        %update Locations
        location(:,:,ii)=cumsum(tangents(:,:,ii)*length_link,2); 
        
        %Compute Distances
        distances(ii,pp)=sum((location(:,end,ii)-location(:,1,ii)).^2);         
    end
    time=toc;
end

%compute theoretical values
theoretical_full=4*length_persist*length_chain-8*length_persist^2*(1-exp(-length_chain/(length_persist*2)));    
theoretical_approx=4*length_persist*length_chain-8*length_persist^2;


%%make non-for loop? more elegant, no need for performance
%%TODO: can you check this? it seems okay now, but id like to have someone
%%check this
for ii=1:P    
    error_chain(ii)=sqrt((mean(distances(:,ii).^2)-mean(distances(:,ii))^2)/N);
end

% plotting figures
close all
if enable_plots
    figure
    errorbar(length_chain,mean(distances,1),error_chain);
    hold on
    plot(length_chain,theoretical_full)
    hold on
    plot(length_chain,theoretical_approx)
    hold off
    title('Length of chain versus the squared end-to-end distance')
    xlabel('Length of chain [nm]')
    ylabel('squared end to end distance [nm^2]')
    legend('Monte Carlo','Theoretical Values','Theoretical Values (using approximation)')
end
%closing statement (for console iterpretability)
fprintf('\n>>> %u configurations completed, Computation finished\n',P)

%% Fitting Function

unknown=linspace(1,500,1000);
% we can determine this using matlab symbolic language, which is probably a
% good thing for bonus points
for ii=1:1000
    Chi_Deriv(ii)=sum(2*(16*unknown(ii)-4*length_chain).*(mean(distances,1)-4*unknown(ii).*length_chain+8*unknown(ii).^2));
end
figure
plot(unknown,Chi_Deriv)
hold on
refline(0,0)


%% Symbolic Language

syms persistence_length_sym
chi=((mean(distances,1)-4*persistence_length_sym.*length_chain+8*persistence_length_sym.^2));
chi_deriv=(2*(16*persistence_length_sym-4*length_chain(end)).*(mean(mean(distances,1))-4*persistence_length_sym*length_chain(end)+8*persistence_length_sym^2));
chi_deriv_2=diff(chi);
chi_deriv_3=persistence_length_sym^2

figure
fplot(chi_deriv)
xlim([1,500])
ylim([-1,5]*10^10)
kappa=root(chi_deriv_2,persistence_length_sym)
vpa(kappa)


%% VECTORISED PERFORMANCE EVALUATION --------------------------------------

