function [straightened_chunks]=straighten_syll_assoc_chunks(syll_assoc_chunks,varargin)
%% Syntax
%
% [straightened_chunks]=straighten_syll_assoc_chunks(syll_assoc_chunks,varargin)
%
%% Inputs  
% syll_assoc_chunks - the syllable asscoiated chunks in the the usual cell
% array format
%
%
%
%% Computation/Processing     
% 
% This function simply unravels the cell array and makes a matrix out of
% it.  
%
% 
%
%% Outputs  
% straightened_chunks - Matrix containing the chunks.
% The rows are frequencies and the columns represent the different
% chunks.
% 
%
%
%% Assumptions
%
%
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
narg_min=1;

if nargin<narg_min
     error(['The number of inputs should at least be ' narg_min])
end

% packaging the inputs into the inputs structure. This can be useful in
% case you need to store the inputs as meta-data with the output. 
inputs=struct('syll_assoc_chunks',syll_assoc_chunks);

% processing supplementary inputs

% Assigning default values to supplementary inputs
supp_inputs.write_to_disk_q=0; % should the function write a file to disk containing its output  
supp_inputs.disk_write_dir='';


supp_inputs=parse_pv_pairs(supp_inputs,varargin);
%

%% Body of the function
no_lengths=size(syll_assoc_chunks{1,1},2);
chosen_length_index=2;
no_freqs_in_chunks=size(syll_assoc_chunks{1,1}{1,chosen_length_index}{1,1},1);
if no_freqs_in_chunks==0
    error('There are no chunks for the chosen length. Choose a different length');
end
straightened_chunks=zeros(no_freqs_in_chunks,1);
last_chunk_no=1;
for i=1:no_lengths
    no_instances=size(syll_assoc_chunks{1,1}{1,i},1);
    no_chunks_per_instance=size(syll_assoc_chunks{1,1}{1,i}{1,1},2);
    if no_chunks_per_instance==0
        continue
    end
    for j=1:no_instances
        straightened_chunks(:,last_chunk_no:last_chunk_no+no_chunks_per_instance-1)=syll_assoc_chunks{1,1}{1,i}{j,1};
        last_chunk_no=last_chunk_no+no_chunks_per_instance;
    end
end
