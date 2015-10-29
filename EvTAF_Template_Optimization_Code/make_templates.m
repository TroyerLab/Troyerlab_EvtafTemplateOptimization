function []=make_templates(target_syll,pre_syll,post_syll,distractor_sylls,data_dir,birdID,requested_lengths,chunk_indices_for_each_requested_length,varargin)
%% Syntax
% 
% []=make_templates(target_syll,pre_syll,post_syll,distractor_sylls,data_dir,birdID,requested_lengths,chunk_indices_for_each_requested_length,varargin)
%
%% Inputs  
%
% target_syll - the symbol of the target syllable for which you want to
% make the template
%
% pre_syll - pre-syllable (if you want to make the template out of instances of the 
% target syllable from a particular sequence)
%
% post_syll - post-syllable (if you want to make the template out of instances of the 
% target syllable from a particular sequence)
% 
% distractor_sylls - symbols of other syllables from the bird's repertoire
% which will serve as distractors
% 
% data_dir - the directory where syll_assoc_chunks file is stored
%
% birdID - the symbol string denoting the bird's identity
%
% requested_lengths - the lengths in chunks of the instances of the target syllable
% from which you want to make the template**1
%
% chunk_indices_for_each_requested_length - the indices of the chunks from
% which you want to make the template**2
% 
%
%% Computation/Processing     
% This function looks for the syll_assoc_chunks file for the syllable and
% the sequence requested. It loads that and picks out all instances of the
% syllable with the requested # of chunks (**1). It then picks up the chunks 
% with the indices given in the input. All these chunks are assembled in
% target_chunks. The template is created by averaging all these chunks.
% 
%
%
% 
%
%% Outputs  
% This function write two files as its output. There is no output provided
% to the calling function. The first file is a .dat file which has the template 
% in the format required by EvTAF. The second file has the same name with
% _metadata.mat appended to the template name. This file contains the metadata for the
% template. Both files are stored in the same directory where the syll_assoc_chunks file
% is found. 
% 
%
%
%% Assumptions
% The template is written by default in the same direcotry as where the
% syll_assoc_chunks file is to be found. 
%
%
%
% % % Triple percentage sign indicates that the code is part of the code
% template and may be activated if necessary in later versions. 
%% Version and Author Identity Notes  
% 
% Last modified by Anand S. Kulkarni on
%
% previous version:
% next version: 
%% Related procedures and functions 
% One needs to run make_syll_assoc_chunks before running this function. 
%
%
%
%% Detailed notes
% **1 - due to variation in duration of syllables, different instances of the syllable
% contain different number of chunks in them. You can pick instances with
% specific lengths to make the template. If you enter 'modal',this function 
% picks the mode of this distribution (of the number of chunks) and works with the instances of
% the target syllable containing the modal number of chunks. this eliminates 
% (somewhat, atleast) the bother about alignment and allows you to average
% over the variation in a given portion of the syllable. 
%
% **2 - if one of your requested length is 12, you may want to make the
% template by using chunks 3 and 4. That is what you would enter here
%
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=6;

prob_path=pwd;

in_message1='Please enter the target syllable';
in_message2='Please enter the symbol of the pre syllable. Hit enter to leave it blank';
in_message3='Please enter the symbol of the post syllable. Hit enter to leave it blank';
in_message4='Please enter the birdID';
in_message5='Please enter the length(s)(in # of chunks) of the target syllable from which you want to make the template. For using the modal length, enter modal. e.g. [5,6,8] or modal';
in_message6='Please enter the index/indices of the chunk that you want to use for making the template. For using the all chunks in a given syllable length, enter all. e.g. {[1,2],''all'',[3,5]}\nNote that this is cell array input';
in_message7='Please select the directory where syllable asscouated chunk files are located';
in_message8='Please enter the symbols of the distractor chunks. e.g. bcdefgh';
if nargin<narg_min
     target_syll=input([in_message1 '\n-->  '],'s');
     pre_syll=input([in_message2 '\n-->  '],'s');
     post_syll=input([in_message3 '\n-->  '],'s');
     distractor_sylls=input([in_message8 '\n-->  '],'s');
     data_dir=uigetdir(prob_path,in_message7);
     birdID=input([in_message4 '\n-->  '],'s');
     requested_lengths=input([in_message5 '\n-->  '],'s');
     chunk_indices_for_each_requested_length=input([in_message6 '\n-->  ']);
end

% packaging the inputs into the inputs structure. This can be useful in
% case you need to store the inputs as meta-data with the output. 
inputs=struct('target_syll',target_syll,'pre_syll',pre_syll,'post_syll',post_syll,'birdID',birdID,'requested_lengths',requested_lengths,'chunk_indices_for_each_requested_length',chunk_indices_for_each_requested_length);

% processing supplementary inputs

% Assigning default values to supplementary inputs
supp_inputs.disk_write_dir=data_dir;
supp_inputs.write_to_disk=1; % should the function write a mat (or for that matter any other data file) file to disk of its output  


supp_inputs=parse_pv_pairs(supp_inputs,varargin);

% Checking if output directory needs to specified and if they have been specified 

if supp_inputs.write_to_disk
    if ~exist(supp_inputs.disk_write_dir,'dir')
        supp_inputs.disk_write_dir=uigetdir('Please select the directory where to store the output files. Hit cancel if you don''t want the function to write a mat file');
        if supp_inputs.disk_write_dir==0
            supp_inputs.write_to_disk=0;
        end
    end
end


% putting file separators at the end of all input paths
if ~strcmpi(data_dir(end),filesep)
    data_dir=[data_dir,filesep];
    supp_inputs.disk_write_dir=[supp_inputs.disk_write_dir,filesep];
end

%

%% Body of the function

% loading up the syll assoc chunks file
seq_str=[lower(pre_syll) upper(target_syll) lower(post_syll)];
syll_assoc_chunks_fullfile=[data_dir 'syll_assoc_chunks_syll_' upper(target_syll) '_seq_' seq_str '.mat'];
if ~exist(syll_assoc_chunks_fullfile,'file')
    error(['A syllable associated chunks file for the ' upper(target_syll) ' does not exist at ' data_dir]);
else
    load(syll_assoc_chunks_fullfile) % loads a variable called syll_assoc_chunks 
end

no_intervals=length(syll_assoc_chunks);

if no_intervals~=1
    error('The number of intervals in the syll_assoc_chunks should be one in order to make a template')
end

no_freqs_in_templ=size(syll_assoc_chunks{1,1}{1,2}{1,1},1);

% obtaining no_instances_for_each_length and lengths
no_lengths=size(syll_assoc_chunks{1,1},2);
no_instances_for_each_length=zeros(1,no_lengths);
lengths=zeros(1,no_lengths);
for i=1:no_lengths
    no_instances_for_each_length(1,i)=size(syll_assoc_chunks{1,1}{1,i},1);
    lengths(1,i)=size(syll_assoc_chunks{1,1}{1,i}{1,1},2);
end

% obtaining indices_for_requested_lengths, requested_lengths (re-obtaining), and no_lengths_requested 
if strcmpi(requested_lengths,'modal')
    [~,max_instances_ind]=max(no_instances_for_each_length);
    modal_length=lengths(max_instances_ind(1));
    requested_lengths=modal_length;
    no_lengths_requested=length(requested_lengths); % this of course, has to be one
    indices_for_requested_lengths=max_instances_ind;
else
    if ischar(requested_lengths)
        requested_lengths=str2num(requested_lengths);
    end
    if isempty(requested_lengths) || ~isnumeric(requested_lengths)
        error('The requested_lengths should be a numeric array')
    end
    no_lengths_requested=length(requested_lengths);
    indices_for_requested_lengths=zeros(1,no_lengths_requested);
    for i=1:no_lengths_requested
        for j=1:no_lengths
            if requested_lengths(1,i)==lengths(1,j)
                indices_for_requested_lengths(1,i)=j;
                break
            end
        end        
    end    
end

% obtaining fully spelled out chunk_indices_for_each_requested_length
if length(requested_lengths)~=length(chunk_indices_for_each_requested_length)
    error('The lengths of no_chunks_in_syll and chunk_index_for_template should be equal');
end

% temp variables for the following operation
temp_indices_for_requested_lengths=indices_for_requested_lengths;
temp_chunk_indices_for_each_requested_length=chunk_indices_for_each_requested_length;

for i=1:no_lengths_requested
   if strcmpi(chunk_indices_for_each_requested_length{1,i},'all')
       no_all_chunks=lengths(indices_for_requested_lengths(1,i));
       temp_indices_for_requested_lengths=[temp_indices_for_requested_lengths(1,1:new_i),repmat(temp_indices_for_requested_lengths(1,new_i),1,no_all_chunks-1),temp_indices_for_requested_lengths(1,new_i+1:end)];
       temp_chunk_indices=cell(1,no_all_chunks);
       for j=1:no_all_chunks
            temp_chunk_indices{1,j}=j;
       end
       temp_chunk_indices_for_each_requested_length=[temp_chunk_indices_for_each_requested_length{1,1:new_i-1},temp_chunk_indices,temp_chunk_indices_for_each_requested_length{1,new_i+1:end}];
       new_i=i+no_all_chunks;
   else
    new_i=i+1;
   end
end

indices_for_requested_lengths=temp_indices_for_requested_lengths;
chunk_indices_for_each_requested_length=temp_chunk_indices_for_each_requested_length;

% over the course of previous code we have obtained: indices_for_requested_lengths
% and  chunk_indices_for_each_requested_length. we will use these going ahead. 
arch_timestamp=datestr(now,'yyyy-mmm-dd HH:MM:SS');
no_lengths_requested=length(indices_for_requested_lengths);
for i=1:no_lengths_requested % equivalent to no_lengths_requested
    no_instances_this_length=no_instances_for_each_length(1,indices_for_requested_lengths(1,i));
    no_target_chunks_per_instance=length(chunk_indices_for_each_requested_length{1,i});       
    target_chunks_for_template=zeros(no_freqs_in_templ,no_instances_this_length*no_target_chunks_per_instance);
       for k=1:no_instances_this_length
         target_chunks_for_template(:,((k-1)*no_target_chunks_per_instance)+1:k*no_target_chunks_per_instance)=syll_assoc_chunks{1,1}{1,indices_for_requested_lengths(1,i)}{k,1}(:,chunk_indices_for_each_requested_length{1,i});
       end
   template=mean(target_chunks_for_template,2);
   template=template-min(template);
   template=template./max(template);

   if supp_inputs.write_to_disk==1
       seq_str=[lower(pre_syll) upper(target_syll) lower(post_syll)];
       template_name=['template_syll_' upper(target_syll) '_seq_' seq_str ...
                      '_chunks_'  num2str(chunk_indices_for_each_requested_length{1,i})...
                      '_outof_' num2str(lengths(indices_for_requested_lengths(1,i)))];
       template_name=strrep(template_name,'  ','_'); % replacing double space with an underscore 
       fid=fopen([data_dir template_name '.dat'],'w');
        for j=1:size(template,1)
            fprintf(fid,'%.5e\n',template(j,1));
        end
        fclose(fid);
        template_metadata.target_syll=target_syll;
        template_metadata.pre_syll=pre_syll;
        template_metadata.post_syll=post_syll;
        template_metadata.distractor_sylls=distractor_sylls;
        template_metadata.file_list=arch_inputs.file_list;
        template_metadata.template=template;
        template_metadata.length_in_chunks_of_target_instances=lengths(indices_for_requested_lengths(1,i));
        template_metadata.target_chunk_indices=chunk_indices_for_each_requested_length{1,i};
        template_metadata.no_instances_used=no_instances_for_each_length(indices_for_requested_lengths(1,i));
        template_metadata.birdID=birdID;
        template_metadata.is_template_synthetic=0;
        template_metadata.synthesis_details=[];
        
        template_metadata.name_of_struct='template_metadata';
        check_struct_against_blank(template_metadata,@create_blank_datastruct_template_metadata);
        
        matfile=[template_name '_metadata.mat'];
        matfullfile=[supp_inputs.disk_write_dir matfile];
        save(matfullfile,'template_metadata','arch_timestamp');
   end    
end


