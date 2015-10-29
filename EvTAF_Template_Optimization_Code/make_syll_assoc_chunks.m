function [syll_assoc_chunks]=make_syll_assoc_chunks(batchfile,batchpath,target_syll,pre_syll,post_syll,varargin)
%% Syntax
%
% [syll_assoc_chunks]=make_syll_assoc_chunks_file(batchpath,batchfile,target_syll,pre_syll,post_syll,varargin)
%
%% Inputs  
%
% batchpath -  location where the batch file is stored
% 
% batchfile -  name of the batch file
% 
% target_syll  - syllable for which you want to make a syllable associated
% chunks file
% 
% pre_syll - pre-syllable (if you want to make chunks out of instances of the 
% target syllable from a particular sequence)
%
% post_syll - post-syllable (if you want to make chunks out of instances of the 
% target syllable from a particular sequence)
% 
%% Computation/Processing     
%  This function will use the *.notmat and *.simdata files for the files listed
%  in the batchfile and find out all instances of the target_syll. For each
% instance, it will pickup all associated chunks. It will sort instances by
% the number of chunks in them and write them in syll_assoc_chunks. 
%
%
% 
%
%% Outputs  
% It writes a syll_assoc_chunks file containing the following cell array. 
% syll_assoc_chunks - is a cell array. it contains within itself a cell array,
% where each cell contains all instances with the same number
% of chunks.   
% 
%
%
%% Assumptions
% Assumes that nomat files and simdata files exist in the location as the batchfile.  
%
%
% % % Triple percentage sign indicates that the code is part of the code
% template and may be activated if necessary in later versions. 
%
%% Version and Author Identity Notes  
% 
% Last modified by Anand S Kulkarni 
%  
%% Related procedures and functions 
% 
%
%
%
%% Detailed notes
% This function is capable of getting chunks associated with multiple intervals. 
% If the target_syll is entered as ab, it will give you all the asscoiated
% chunks for three intervals a, gap between a and b and b. 
%
%
%
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=5;

prob_path=pwd;

in_message1='Please select the batch file with the list of files you want to use';
in_message2='Please enter the symbol of the target syllable';
in_message3='Please enter the symbol of the pre syllable. Hit enter to leave it blank';
in_message4='Please enter the symbol of the post syllable. Hit enter to leave it blank';
if nargin<narg_min
     [batchfile,batchpath]=uigetfile([prob_path filesep 'batch*'],in_message1);    
     target_syll=input([in_message2 '\n-->  '],'s'); 
     pre_syll=input([in_message3 '\n-->  '],'s'); 
     post_syll=input([in_message4 '\n-->  '],'s'); 
end

% packaging the inputs into the inputs structure. This can be useful in
% case you need to store the inputs as meta-data with the output. 
inputs=struct('batchfile',batchfile,'batchpath',batchpath,'target_syll',target_syll,'pre_syll',pre_syll,'post_syll',post_syll);


% processing supplementary inputs

% Assigning default values to supplementary inputs
supp_inputs.disk_write_dir=batchpath;
supp_inputs.write_to_disk_q=1; % should the function write a mat file to disk containing its output  

supp_inputs=parse_pv_pairs(supp_inputs,varargin);

% Checking if output directories need to specified and if they have been specified 

if supp_inputs.write_to_disk_q
    if ~exist(supp_inputs.disk_write_dir,'dir')
        supp_inputs.disk_write_dir=uigetdir('Please select the directory where to store the output mat file. Hit cancel if you don''t want the function to write a mat file');
        if supp_inputs.disk_write_dir==0
            supp_inputs.write_to_disk_q=0;
        end
    end
end


% putting file separators at the end of all input paths
if ~strcmpi(batchpath(end),filesep)
    batchpath=[batchpath,filesep];
    supp_inputs.disk_write_dir=[supp_inputs.disk_write_dir,filesep];
end



%
%% Body of the function

fid=fopen([batchpath batchfile],'r');
file=fgetl(fid);
chunks=cell(0);
no_chunks_vec=[];
file_list=cell(0);

no_intervals=2*length(target_syll)-1;
srch_string=[pre_syll target_syll post_syll];

while ischar(file)
    file
    [~,name,ext]=fileparts(file);    
    notmat_file=[batchpath name ext '.not.mat'];
    simdata_file=[batchpath name '_simdata.mat'];
    load(notmat_file)
    load(simdata_file)
    file_list=[file_list;file];
    
   inds=strfind(labels,srch_string);
   no_isntances_in_file=length(inds);
   temp_chunks=cell(length(inds),no_intervals);
   temp_no_chunks_vec=zeros(length(inds),no_intervals);

   for j=1:no_isntances_in_file

       for l=1:no_intervals

           if mod(l,2)==0
               interval_ind=[inds(j)+l/2-1,inds(j)+l/2]+length(pre_syll);
           else
               interval_ind=inds(j)+((l+1)/2)-1+length(pre_syll);
           end

           assoc_chunks=[];

           for k=1:length(simdata.chunk_assocs)
               if (simdata.chunk_assocs{k})==interval_ind 
                   assoc_chunks=[assoc_chunks,k];
               end
           end 

          temp_chunks{j,l}=simdata.chunks(:,assoc_chunks);
          temp_no_chunks_vec(j,l)=length(assoc_chunks);         
       end
   end
   chunks=[chunks;temp_chunks];
   no_chunks_vec=[no_chunks_vec;temp_no_chunks_vec];
   file=fgetl(fid) ;
end

fclose(fid);

no_instances=size(chunks,1);
unique_lengths=cell(1,no_intervals);
syll_assoc_chunks=cell(1,no_intervals);

for j=1:no_intervals

    [lengths_vec,~]=count_unique(no_chunks_vec(:,j));
    unique_lengths{1,j}=lengths_vec;    

    syll_assoc_chunks{1,j}=cell(1,length(unique_lengths{1,j}));

    for k=1:no_instances 
        instance_chunks=chunks{k,j};
        no_chunks_in_instance=size(instance_chunks,2);

        len_index=find(unique_lengths{1,j}==no_chunks_in_instance);
        syll_assoc_chunks{1,j}{1,len_index}=[syll_assoc_chunks{1,j}{1,len_index};{instance_chunks}];
    end
end


%
%% Processing outputs and ending stuff
arch_timestamp=datestr(now,'yyyy-mmm-dd HH:MM:SS');
inputs.file_list=file_list;
arch_inputs=inputs;
arch_supp_inputs=supp_inputs;

if supp_inputs.write_to_disk_q==1
    seq_str=[lower(pre_syll) upper(target_syll) lower(post_syll)];
    matfile=['syll_assoc_chunks_syll_' upper(target_syll) '_seq_' seq_str '.mat'];
    matfullfile=[supp_inputs.disk_write_dir matfile];
    save(matfullfile,'syll_assoc_chunks','arch_inputs','arch_supp_inputs','arch_timestamp');
end


% % % if supp_inputs.write_fig_q==1
% % %     figfile='function_plot.fig';
% % %     figfullfile={[supp_inputs.disk_write_dir filesep figfile]};
% % %     saveas(fig1,figfullfile);
% % % end

% removing the stop that was put for easier debugging
dbclear if error



