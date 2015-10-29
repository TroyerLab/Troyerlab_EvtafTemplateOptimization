function [results]=optimize_template(template_file,template_path,target_input_method,...
                  target_file,target_path,distractor_input_method,distractor_file,...
                  distractor_path,distractor_factor,varargin)
%% Syntax
% [results]=optimize_template(template_file,template_path,target_input_method,...
%                   target_file,target_path,distractor_input_method,distractor_file,...
%                   distractor_path,distractor_factor,varargin)
%              
%% Inputs and Supplemental Inputs  
%
% template_file - name of the template file you want to optimize
%
% template_path - location of the template file
%
% target_input_method - where should the function find the target chunks. two
% options. 'template_metadata': load the syll associated chunks for the length and 
% the index listed in the template metadata. 'file_input' :load 
% the target chunks from a given mat file.
%
% if input_method = 'file_input'
% target_file - name of the file containing target chunks
% 
% target_path - location of the file containing target chunks 
%
% distractor_input_method -where should the function find the distractor chunks. two
% options. 'template_metadata': load the syll associated chunks for all the 
% listed distractor syllables in the template metadata. 'file_input' :load 
% the distractor chunks from a given mat file.
%
% if input_method = 'file_input'
% distractor_file - name of the file containing distractor chunks
% 
% distractor_path - location of the file containing distractor chunks
%
% distractor_factor - the factor by which the probability density function
% of the distractor distances should be multiplied **1%
%
% supp_inputs.sigma_init - the initial sigma of the Gaussians used to
% construct the continuous distributions.default = 0.5;
% 
% supp_inputs.collation_function=function used to collate individual
% gaussians and construct the continuous density distributions. default = @mean.
% An alternative can be @sum. 
%
% supp_inputs.target_factor = the factor by which the probability density function
% of the target distances is multipled. similar to distractor_factor.**2. default=1.   
%
% supp_inputs.template_step_mag - The size of the step for each iteration of
% the optimization.
%
%% Computation   
% 
% Modifies the template in a manner dictated by the gradient descent algorithm
% to achieve minimum overlap (error) between the distance distributions of
% target and distractor chunks. 
%
%
% 
%
%% Outputs  
%
% optim_archive.error_vec= vector containing magnitude of error through all
% the steps of the optimization
%
% optim_archive.threshold_vec=vector containing value of the threshold through all
% the steps of the optimization
%
% optim_archive.fne_instances_vec=vector containing # of false negative chunks through all
% the steps of the optimization
%
% optim_archive.fpe_instances_vec=vector containing # of false positive chunks through all
% the steps of the optimization
%
% optim_archive.fne_rate_vec=vector containing % of false negative chunks through all
% the steps of the optimization
%
% optim_archive.fpe_rate_vec=vector containing % of false positive chunks through all
% the steps of the optimization
%
% optim_archive.fin_sigma_vec=vector containing value of the sigma used through all
% the steps of the optimization
%
% optim_archive.template_grad_vec=matrix containing the gradient vector through all
% the steps of the optimization
%
% optim_archive.template_vec=matrix containing the template vector through all
% the steps of the optimization
% 
% optim_archive.template_grad_mag_vec=vector containing magnitude of the gradient of the template
% through all the steps of the optimization
%
%% Assumptions
%
%
%
%
%
%% Version and Author Identity Notes  
% Last modified by The Anand S Kulkarni on 
%
% previous version:
% next version: 
%% Related procedures and functions 
% 
%
%
%
%% Detailed Notes
% 
% %**1 - this effectively determines where the threshold gets placed by
% shifting the crossing point of the two distributions (in the valley). lower 
% values shift the crossing point to the right and higher values shift it
% to the left. It also determines the relative contribution (relative to target chunks)
% of the distractors to the gradient and the error calculation. 
% 
%%**2 - the target factor is by default kept at 1. Moving around of threhsold is
% achieved by controlling the distractor_factor itself. Other aspects involving relative 
% contributions of target and distractor chunks are also similarly
% controlled. 
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=9;

prob_path=pwd;

in_message1='Please select the template file which you want to optimize';
in_message2=['Please enter the input method by which target chunks should be'...
            'located. ''template_metadata'' or ''file_input'''];
in_message3='Please select the file containing target chunks';
in_message4=['Please enter the input method by which distractor chunks should be'...
            'located. ''template_metadata'' or ''file_input'''];
in_message5='Please select the file containing distractor chunks';
in_message6='Please enter the the factor by which the probability density function of the distractor distances should be multiplied while determining the threshold';
if nargin<narg_min
   [template_file,template_path]=uigetfile([prob_path filesep '*.dat'],in_message1); 
   
   target_input_method=input([in_message2 '\n-->  '],'s'); 
     if strcmpi(target_input_method,'file_input') 
        [target_file,target_path]=uigetfile([prob_path filesep '*.mat'],in_message3);  
     elseif strcmpi(target_input_method,'template_metadata') 
         target_file='';
         target_path='';
     else
         error('Incorrect input for the variable target_input_method. Likely a typo.')
     end
     
     distractor_input_method=input([in_message4 '\n-->  '],'s'); 
     if strcmpi(distractor_input_method,'file_input') 
        [distractor_file,distractor_path]=uigetfile([prob_path filesep '*.mat'],in_message5);  
     elseif strcmpi(distractor_input_method,'template_metadata') 
         distractor_file='';
         distractor_path='';
     else
         error('Incorrect input for the variable distractor_input_method. Likely a typo.')
     end
     
     distractor_factor=input([in_message6 '\n-->  ']); 
     
end

% packaging the inputs into the inputs structure. This can be useful in
% case you need to store the inputs as meta-data with the output. 
inputs=struct('template_file',template_file,'template_path',template_path...
              ,'target_input_method',target_input_method,'target_file',target_file...
              ,'target_path',target_path,'distractor_input_method',distractor_input_method...
              ,'distractor_file',distractor_file,'distractor_path',distractor_path...
              ,'distractor_factor',distractor_factor);

% processing supplementary inputs

% Assigning default values to supplementary inputs
supp_inputs.sigma_init=0.5;
supp_inputs.collation_function=@mean;
supp_inputs.target_factor=1;
supp_inputs.template_step_mag=5;
supp_inputs.write_to_disk_q=0; % should the function write a file to disk containing its output  
supp_inputs.disk_write_dir='';


supp_inputs=parse_pv_pairs(supp_inputs,varargin);

% putting file separators at the end of all input paths

if ~isempty(template_path)
   if ~strcmpi(template_path(end),filesep)
        template_path=[template_path,filesep];
   end
end
if ~isempty(target_path)
    if ~strcmpi(target_path(end),filesep)
        target_path=[target_path,filesep];
    end
end
if ~isempty(distractor_path)
    if ~strcmpi(distractor_path(end),filesep)
        distractor_path=[distractor_path,filesep];
    end
end

%

%% Body of the function

% loading template, target_chunks, and distractor_chunks

template=load([template_path template_file]);

[target_chunks]=assemble_target_chunks(template_file,template_path,target_input_method,...
                  target_file,target_path);
[distractor_chunks]=assemble_distractor_chunks(template_file,template_path,distractor_input_method,...
                  distractor_file,distractor_path);


% initiating figure and axes
grand_fig=figure;
densities_axes=subplot(1,3,1,'parent',grand_fig);
error_axes=subplot(1,3,2,'parent',grand_fig);
t_grad_mag_axes=subplot(1,3,3,'parent',grand_fig);

% initaitng archive variables
optim_archive.error_vec=[];
optim_archive.threshold_vec=[];
optim_archive.fne_instances_vec=[];
optim_archive.fpe_instances_vec=[];
optim_archive.fne_rate_vec=[];
optim_archive.fpe_rate_vec=[];
optim_archive.fin_sigma_vec=[];
optim_archive.template_grad_vec=[];
optim_archive.template_vec=[];
optim_archive.template_grad_mag_vec=[];

% other initializations
first_entry=1;
iteration_no=0;

while true 
    iteration_no=iteration_no+1
    iter_vec=(1:iteration_no);

    % calculating distances
    dists_target=sqrt(sum((target_chunks-repmat(template,1,size(target_chunks,2))).^2,1));
    dists_distractor=sqrt(sum((distractor_chunks-repmat(template,1,size(distractor_chunks,2))).^2,1));

    % calculating optimal sigma/ verifying optimal sigma
    [sigma,density_target,density_distractors,dists_vec]...
    =calculate_optimal_sigma(dists_target,dists_distractor,'verify_only',1,...
    'sigma_init',supp_inputs.sigma_init);

    % setting the sigma_init of the next iteration as the current optimal
    % sigma
    supp_inputs.sigma_init=sigma;
    sigma_target=sigma;
    sigma_distractor=sigma;
    
    % calculating optimal threshold
    [threshold]=calculate_optimal_threshold(density_target,density_distractors,dists_vec,distractor_factor); 
    
     % multiplying the density distributions with appropriate factors for
    % plotting
    density_target=density_target*supp_inputs.target_factor;
    density_distractors=density_distractors*distractor_factor;
       
    
    % archiving threshold and sigma
    optim_archive.threshold_vec=[optim_archive.threshold_vec,threshold];
    optim_archive.fin_sigma_vec=[optim_archive.fin_sigma_vec,sigma];
    
   
    
    % plotting densities and threshold    
    if first_entry==1 
        densities_in_h=plot(densities_axes,dists_vec,density_target,'g');
        hold(densities_axes,'on')
        densities_out_h=plot(densities_axes,dists_vec,density_distractors,'r');
        th_ylim=[0,max([density_target;density_distractors])];
        th_xlim=[threshold,threshold];
        th_line_h=plot(densities_axes,th_xlim,th_ylim,'color','k','linestyle','--');
        hold(densities_axes,'off')
        title(densities_axes,'Density of Distance Distributions and the Threshold')
        set(densities_in_h,'XDataSource','dists_vec','YdataSource','density_target');
        set(densities_out_h,'XDataSource','dists_vec','YdataSource','density_distractors');
        set(th_line_h,'Xdatasource','th_xlim','Ydatasource','th_ylim');
    else
        % refreshing
        th_ylim=[0,max([density_target;density_distractors])];
        th_xlim=[threshold,threshold];
    end
    

   
    %  error rates and instances: calculation and archival
    optim_archive.fne_instances_vec=[optim_archive.fne_instances_vec,length(find(dists_target>=threshold))];
    optim_archive.fpe_instances_vec=[ optim_archive.fpe_instances_vec,length(find(dists_distractor<=threshold))];
    optim_archive.fne_rate_vec=[optim_archive.fne_rate_vec,length(find(dists_target>=threshold))/length(dists_target)];
    optim_archive.fpe_rate_vec=[optim_archive.fpe_rate_vec,length(find(dists_distractor<=threshold))/length(dists_distractor)];

    
    % calculating total error and archiving it
    fn_error=supp_inputs.collation_function(1-cdf('norm',threshold.*ones(size(dists_target)),dists_target,sigma_target*ones(size(dists_target))));
    fp_error=supp_inputs.collation_function(cdf('norm',threshold.*ones(size(dists_distractor)),dists_distractor,sigma_distractor*ones(size(dists_distractor))));
    total_error=supp_inputs.target_factor*fn_error+distractor_factor*fp_error;
    optim_archive.error_vec=[optim_archive.error_vec,total_error];
    
    % plotting total error
    if first_entry
        err_line_h=plot(error_axes,optim_archive.error_vec,'*-'); 
        title(error_axes,'Total Error')
        set(err_line_h,'ydatasource','optim_archive.error_vec','xdatasource','iter_vec');
    end   

    
    % template gradient calculation
    template_grad_fne=supp_inputs.collation_function((repmat(template,1,size(target_chunks,2))-target_chunks) .*repmat(pdf('norm',threshold*ones(size(dists_target)),...
        dists_target,sigma_target*ones(size(dists_target)))./dists_target,size(target_chunks,1),1),2);% fne=false negative error
    
    template_grad_fpe=supp_inputs.collation_function((repmat(template,1,size(distractor_chunks,2))-distractor_chunks) .*repmat(pdf('norm',threshold*ones(size(dists_distractor)),...
         dists_distractor,sigma_distractor*ones(size(dists_distractor)))./dists_distractor,size(distractor_chunks,1),1),2);% fpe=false positive error
     
    template_grad=supp_inputs.target_factor*template_grad_fne-distractor_factor*template_grad_fpe;
    template_grad_mag=sqrt(sum(template_grad.^2));
    
    % archiving template gradient and its magnitude
    optim_archive.template_grad_vec=[optim_archive.template_grad_vec,template_grad];  
    optim_archive.template_grad_mag_vec=[optim_archive.template_grad_mag_vec,template_grad_mag];

    % plotting the magnitude of the gradient of the template
    if first_entry       
       template_grad_mag_h=plot(t_grad_mag_axes,optim_archive.error_vec,'*-'); 
       title(t_grad_mag_axes,'Template Gradient Magnitude')
       set(template_grad_mag_h,'ydatasource','optim_archive.template_grad_mag_vec','xdatasource','iter_vec');
    end    
    

    % minimum reached check, refreshing all plots, and updating the template  
    if first_entry
        first_entry=0;    
        template=template-supp_inputs.template_step_mag*template_grad;
        optim_archive.template_vec=[optim_archive.template_vec,template];
    else
        refreshdata(grand_fig,'caller')
        drawnow  
        if optim_archive.error_vec(end-1)/optim_archive.error_vec(end)<1.0001 && optim_archive.error_vec(end-1)/optim_archive.error_vec(end)>1.0
                break
        elseif optim_archive.error_vec(end-1)/optim_archive.error_vec(end)<1.0
            supp_inputs.template_step_mag=supp_inputs.template_step_mag/2;
            template=template-supp_inputs.template_step_mag*template_grad;
            optim_archive.template_vec=[optim_archive.template_vec,template];
        else % get a new template
            template=template-supp_inputs.template_step_mag*template_grad;
            optim_archive.template_vec=[optim_archive.template_vec,template];
        end       
    end
    
end

%% Processing outputs 
arch_timestamp=datestr(now,'yyyy-mmm-dd HH:MM:SS');
arch_inputs=inputs;
arch_supp_inputs=supp_inputs;

results.arch_timestamp=arch_timestamp;
results.arch_inputs=arch_inputs;
results.arch_supp_inputs=arch_supp_inputs;
results.optim_archive=optim_archive;


% alternaitve code for the last portion (NOT COMPREHENSIVELY TESTED)

% % % if optim_archive.error_vec(end-1)/optim_archive.error_vec(end)<1.0001 && optim_archive.error_vec(end-1)/optim_archive.error_vec(end)>1.0
% % %     if optim_archive.error_vec(end)==min(optim_archive.error_vec)
% % %         break
% % %     else
% % %         [~,min_ind]=min(optim_archive.error_vec);
% % %         template=optim_archive.template_vec(:,min_ind)-step_fac_t*t_grad_vec(:,min_ind);
% % %         optim_archive.template_vec=[optim_archive.template_vec,template];
% % %     end
% % % elseif optim_archive.error_vec(end-1)/optim_archive.error_vec(end)<1.0
% % %         supp_inputs.template_step_mag=supp_inputs.template_step_mag/2;
% % %         template=template-supp_inputs.template_step_mag*template_grad;
% % %         optim_archive.template_vec=[optim_archive.template_vec,template];
% % % else % get a new template
% % %     template=template-supp_inputs.template_step_mag*template_grad;
% % %     optim_archive.template_vec=[optim_archive.template_vec,template];
% % % end
        







    

