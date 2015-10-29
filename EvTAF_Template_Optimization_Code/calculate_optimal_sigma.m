function [sigma,density_target,density_distractors,dists_vec]=calculate_optimal_sigma(dists_target,dists_distractor,varargin)
%% Syntax
%
% [sigma,density_target,density_distractors,dists_vec]=calculate_optimal_sigma(dists_target,dists_distractor,varargin)
%
%% Inputs  
%
% dists_target - vector containing the distances for target chunks
%
% dists_distractor - vector containing the distances for distractor chunks
% 
% supp_inputs.sigma_init - the initial sigma of the Gaussians used to
% construct the continuous distributions.default = 0.5;
%
% supp_inputs.sigma_step - the quantum of increment or decrement in sigma. 
% can used while finding optimal sigma. default = 0.05;   
%
% supp_inputs.collation_function= function used to collate individual
% gaussians and construct the continuous density distributions. default = @mean.
% An alternative can be @sum. 
%
% supp_inputs.no_bins - the number of bins into which the distance axis
% needs to be divided (between min and max distances). default=1000. 
%
% supp_inputs.verify_only -  if this is 1, the function will simply check
% if, for the given value of sigma_init, the density distributions are
% monotonous (to the right of the highest peak for target density distribution
% and to the left of the highest peak for the distractor density
% distribution). If they are, it will simply return the sigma init as the output sigma.  
% default = 0. 
%
% supp_inputs.case_verify_failed- this input controls what happens if the
% verification fails (either one or both distributions are not monotonous).
% if this is set to 'throw_error', the function will throw an error. If it
% is set to anything else, the function will find the optimal sigma and
% output it. 
% default = 'throw_error' 
% 
% 
%% Computation/Processing     
% 
% This function calculates the optimal sigma (standard deviation for Gaussians)
% used for constructing the continous distributions for target and distractor distances.
% It can also verify that a given sigma, supplied to it, is optimal or not.
% 
% An optimal sigma is one which gives monotonous distributions for both target 
% and distractor chunks. The target distributions is required to be monotonous 
% to the right of its highest peak and, likewise, the distractor
% distribution is required to be monotonous to he left of its highest peak.
% 
%
% 
%
%% Outputs  
% sigma - optimal sigma computed or verified by the function 
% 
% density_target -  the probability density distribution for target distances 
% 
% density_distractors -the probability density distribution for distractor distances 
%
% dists_vec - the domain on which the above two distributions is
% calculated. this is the discretized distance axis. 
% 
%
%
%% Assumptions
%
%
%
% % % Triple percentage sign indicates that the code is part of the code
% template and may be activated if necessary in later versions. 
%% Version and Author Identity Notes  
% 
% Last modified by Anand S Kulkarni on 
% 
% previous version:
% next version: 
%% Related procedures and functions 
% 
%
%
%
%% Detailed notes
%
%
%
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=3;

if nargin<narg_min
     error(['The number of inputs should at least be ' narg_min])
end

% processing supplementary inputs

% Assigning default values to supplementary inputs
supp_inputs.sigma_init=0.5;
supp_inputs.sigma_step=0.05;
supp_inputs.no_bins=1000;
supp_inputs.collation_function=@mean;
supp_inputs.verify_only=0;
supp_inputs.case_verify_failed='throw_error';
supp_inputs.write_to_disk_q=0; % should the function write a file to disk containing its output  
supp_inputs.disk_write_dir='';


supp_inputs=parse_pv_pairs(supp_inputs,varargin);

%

%% Body of the function

% discretizing the distance axis
dist_min=min(min(dists_target),min(dists_distractor));
dist_max=max(max(dists_target),max(dists_distractor));
dist_incre=(dist_max-dist_min)/supp_inputs.no_bins;
dists_vec=(dist_min:dist_incre:dist_max)';

% standard deviations of gaussians 
sigma=supp_inputs.sigma_init;
sigma_target=sigma;
sigma_distractors=sigma;    

% determining appropriate value of sigma
first_entry=1;
direction='up';
exit_condition_met=0;
while ~exit_condition_met
 
    % incrementing or decrementing sigma for a possible next iteration
    if ~first_entry
        if strcmpi(direction,'down')
            sigma_prev=sigma;
            sigma=sigma-supp_inputs.sigma_step;
        else
           sigma=sigma+supp_inputs.sigma_step;
        end    
        sigma_target=sigma;
        sigma_distractors=sigma;
    end

    % calculating densities
    density_target=supp_inputs.collation_function...
                (pdf('norm',repmat(dists_vec,1,size(dists_target,2)),repmat...
                (dists_target,size(dists_vec,1),1),sigma_target*ones(size(dists_vec,1),size(dists_target,2))),2);
            
    density_distractors=supp_inputs.collation_function...
                  (pdf('norm',repmat(dists_vec,1,size(dists_distractor,2)),repmat...
                  (dists_distractor,size(dists_vec,1),1),sigma_distractors*ones(size(dists_vec,1),size(dists_distractor,2))),2);

    % calculating # of inflections to the right of the highest one 

    % determining the highest inflection
    % the plus one indicates that if the fifth point was <0, it means the sixth point has an inflection 
    inflections_target_inds=find(exp(diff(log(diff(density_target))))<0)+1;
    inflections_distractors_inds=find(exp(diff(log(diff(density_distractors))))<0)+1;
    [~,max_infection_target_ind]=max(density_target(inflections_target_inds));
    [~,max_infection_distractors_ind]=max(density_distractors(inflections_distractors_inds));
    
    

    if max_infection_target_ind==length(inflections_target_inds) % if the inflection with the max value is the last one, we are good
        target_go=0;
    else
        target_go=1;
    end
    if max_infection_distractors_ind==1 % if the inflection with the max value is the first one, we are good
        distractors_go=0;
    else
        distractors_go=1;
    end
    
    if first_entry
         if ~target_go&&~distractors_go 
             if supp_inputs.verify_only 
                 sigma_verified=1;
                 return
             else
                 direction='down';
             end
         else
             if supp_inputs.verify_only 
                 if strcmpi(supp_inputs.case_verify_failed,'throw_error')
                     error('The sigma supplied is not sufficient for constructing smooth, monotonous density distributions')                    
                 end             
             end       
             
         end
         first_entry=0;
     end

    
    if strcmpi(direction,'down')
        exit_condition_met=target_go||distractors_go;
    else
       exit_condition_met=~target_go&&~distractors_go;
    end       
end

if strcmpi(direction,'down')
    sigma=sigma_prev;
    density_target=supp_inputs.collation_function...
            (pdf('norm',repmat(dists_vec,1,size(dists_target,2)),repmat...
            (dists_target,size(dists_vec,1),1),sigma_target*ones(size(dists_vec,1),size(dists_target,2))),2);
            
    density_distractors=supp_inputs.collation_function...
                  (pdf('norm',repmat(dists_vec,1,size(dists_distractor,2)),repmat...
                  (dists_distractor,size(dists_vec,1),1),sigma_distractors*ones(size(dists_vec,1),size(dists_distractor,2))),2);
end