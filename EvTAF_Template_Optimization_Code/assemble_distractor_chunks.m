function [distractor_chunks]=assemble_distractor_chunks(template_file,template_path,...
                   input_method,distract_file,distract_path,varargin)
%% Syntax
%
% [distractor_chunks]=assemble_distractor_chunks(template_file,template_path,...
%                   input_method,distract_file,distract_path,varargin)
%% Inputs  
%
% template_file - name of the template file
% 
% template_path - location of the template file
%
% input_method - where should the function find the distractor chunks. two
% options. 'template_metadata': load the syll associated chunks for all the 
% listed distractor syllables in the template metadata. 'file_input' :load 
% the distractor chunks from a given mat file.
%
% if input_method = 'file_input'
% distract_file - name of the file containing distractor chunks
% 
% distract_path - location of the file containing distractor chunks
%
% supp_inputs.missing_distract_files -  what to do in case the syll assoc
% files for some of the distractor syllables give an error while loading.
% two options: 'throw_error' - this is the default value. the code will throw an
% error. 'ignore' - the code will ignore the error encountered while
% loading the file and will continue to the next distractor syllable in the
% list. 
%
% supp_inputs.distract_verify - the file containing the distractor chunks 
% may contain some metadata (depends on if the user put it there while making that file). 
% this metadata allows the function to verify that the distractor chunks are 
% indeed meant for the template being referred to at the beginning of the function. 
% if this field is set to 1 (the default value), the function will try to
% perform that verification. if it is set to 0, the function will not
% verify. 
%
%% Computation/Processing     
% assembles distractor chunks into variable called distractor_chunks
%
%
% 
%
%% Outputs  
%  distractor_chunks -  matrix containing the distractor chunks. the rows
%  represent frequencies and the columns are the different chunks. 
% 
%
%
%% Assumptions
%
%  The function assumes that the syll_assoc_chunks are located at the same
%  path as the template
%
%
% % % Triple percentage sign indicates that the code is part of the code
% template and may be activated if necessary in later versions. 
%% Version and Author Identity Notes  
% 
% Last modified by Anand Kulkarni on 
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
%
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=5;

prob_path=pwd;

in_message1='Please select the template file';
in_message2=['Please enter the input method by which distractor chunks should be'...
            'located. ''template_metadata'' or ''file_input'''];
in_message3='Please select the file containing distractor chunks';
if nargin<narg_min
     [template_file,template_path]=uigetfile([prob_path filesep '*.dat'],in_message1);   
     input_method=input([in_message2 '\n-->  '],'s'); 
     if strcmpi(input_method,'file_input') 
        [distract_file,distract_path]=uigetfile([prob_path filesep '*.mat'],in_message3);  
     elseif strcmpi(input_method,'template_metadata') 
         distract_file='';
         distract_path='';
     else
         error('Incorrect input for the variable input_method. Likely a typo.')
     end
end

% packaging the inputs into the inputs structure. This can be useful in
% case you need to store the inputs as meta-data with the output. 
inputs=struct('template_file',template_file,'template_path',template_path,...
     'input_method',input_method,'distract_file',distract_file,'distract_path',distract_path);

% processing supplementary inputs

% Assigning default values to supplementary inputs
supp_inputs.missing_distract_files='throw_error';
supp_inputs.distract_verify=1;
supp_inputs.write_to_disk_q=0; % should the function write a file to disk containing its output  
supp_inputs.disk_write_dir='';

supp_inputs=parse_pv_pairs(supp_inputs,varargin);

% Checking if output directories need to specified and if they have been specified 

if supp_inputs.write_to_disk_q
    if ~exist(supp_inputs.disk_write_dir,'dir')
        supp_inputs.disk_write_dir=uigetdir(['Please select the directory where'...
                                   'to store the output mat file. Hit cancel if'... 
                                   'you don''t want the function to write a mat file']);
        if supp_inputs.disk_write_dir==0
            supp_inputs.write_to_disk_q=0;
        end
    end
end

% putting file separators at the end of all input paths
if ~strcmpi(template_path(end),filesep)
    template_path=[template_path,filesep];
end

if ~isempty(distract_path)
    if ~strcmpi(distract_path(end),filesep)
        distract_path=[distract_path,filesep];
    end
end

if supp_inputs.write_to_disk_q
    if ~strcmpi(supp_inputs.disk_write_dir(end),filesep)
        supp_inputs.disk_write_dir=[supp_inputs.disk_write_dir,filesep];
    end
end

%% Body of the function

% loading template metadata
[~,metadata_filename,~]=fileparts(template_file);
metadata_fullfile=[template_path metadata_filename '_metadata.mat'];
load(metadata_fullfile) % loads a variable called template_metadata
no_freqs_in_templ=size(template_metadata.template,1);

if strcmpi(input_method,'template_metadata')
    
    distractor_sylls=template_metadata.distractor_sylls;
    no_distract_sylls=size(distractor_sylls,2);
    distractor_chunks=zeros(no_freqs_in_templ,1);
    last_chunk_no=1;
    % loading distractor_chunks
    for i=1:no_distract_sylls
        distr_syll=distractor_sylls(1,i);
        syll_assoc_chunks_fullfile=[template_path 'syll_assoc_chunks_syll_'...
                                upper(distr_syll) '_seq_' upper(distr_syll) '.mat'];
        try 
            load(syll_assoc_chunks_fullfile); % loads a bunch of variables. syll_assoc_chunks is the one we need
        catch exp
            if strcmpi(exp.identifier,'MATLAB:load:couldNotReadFile')
                if strcmpi(supp_inputs.missing_distract_files,'ignore')                    
                    continue
                else
                    rethrow(exp);
                end                
            else
               rethrow(exp);
            end
        end
        
        no_freqs_in_syll_assoc_chunks=size(syll_assoc_chunks{1,1}{1,2}{1,1},1);
        if no_freqs_in_syll_assoc_chunks~=no_freqs_in_templ
             error('The number of freqs in the loaded chunks and in the template do not match')
        end
        
        straightened_chunks=straighten_syll_assoc_chunks(syll_assoc_chunks);
        no_st_chunks=size(straightened_chunks,2);
        distractor_chunks(:,last_chunk_no:last_chunk_no+no_st_chunks-1)=straightened_chunks;
        last_chunk_no=last_chunk_no+no_st_chunks;
    end
    
elseif strcmpi(input_method,'file_input')
    
    if exist([distract_path distract_file],'file')
        load([distract_path distract_file]) % loads a variable called distractor_chunks
    else
       [distract_file,distract_path]=uigetfile([prob_path filesep '*.dat'],in_message3);  
        load([distract_path distract_file]) % loads a variable called distractor_chunks
    end
    % checking if distractor_chunks have been loaded correclty
    if ~exist('distractor_chunks','var')
        error('The loaded file does not contain varaible called distractor chunks')
    else
        if size(distractor_chunks,1)~=no_freqs_in_templ
            error('The number of freqs in the loaded chunks and in the template do not match')
        else
            if supp_inputs.distract_verify
               % palce holder for other checks to be written later. this will be done one 
               % the exact structure of the distractor_chunks meta-data is decided 
            end            
        end
    end
    
else
  error('Incorrect input for the variable input_method. Likely a typo.')    
end
