function []=create_simdata_file(batchname,batchpath,chanspec,pretime,fs,use_audio_metadata,templname,templpath,varargin)
%% Syntax
%
% []=create_simdata_file(batchname,batchpath,chanspec,pretime,fs,use_audio_metadata,templname,templpath,varargin)
%
%% Inputs  
%
% batchname - name of the batch file containing the cbin files
% 
% batchpath - location of the batch file
%
% chanspec -  the spcification of the channel to be read from the cbin
% files
% 
% pretime - the pre-time (in seconds) that you want to apply these audio files as
% default. it will be used if use_audio_metadata is not available. 
%
% fs - the sampling frequency (in Hz)  that you want to apply these audio files as
% default. it will be used if use_audio_metadata is not available. 
% 
% use_audio_metadata - do you want to use audio metadata like the template used during
% recording, the tmp and the rec files generated during recording
% 
% if use_audio_metadata=1
% templname - name of the template used during recording
% 
% templpath - location of the template
% 
% supp_inputs.window_min - the number of offsets on either side of the offset
% with the minimum difference between simulated and real distances. default=20; 
%
% supp_inputs.overwrite -  if a given audio file already has a simdata file, should the
% function overwrite it. default='n';
%
% supp_inputs.no_dists_for_compar - the number of distances used for evaluating 
% a given offset setting. default=500;
% 
%% Computation/Processing     
% 
% If use_audio_metadata =1 
% It creates the *_simdata.mat file for each cbin file. this file contains
% information on the first sample point that evtaf used for chunking the
% incoming data stream. it does so by comparing the simulated distance
% values to real ones (obtained from the tmp files). in order to find out
% the exact point, it does this comparison at the sample specified by fix(pretime*fs) 
%  and samples offset from that point by a range of values. the offset at   
% fix(pretime*fs) is of course zero. it finds out the offset that gives you
% the minimum difference between simulated and real distance values. it
% uses the sample point specified by that offset as the first point of the
% first chunks in the file. 
% 
% If use_audio_metadata =0 
% it simply assumes the first point of the first chunk to be at fix(pretime*fs) 
% 
%
%% Outputs  
% 
% the function saves a file *_simdata.mat for each cbin file it reads at
% the same location as the batch file. it contains a struct called simdata 
% with the following fields:
%
% simdata.first_chunk_first_point -  the first point of the first chunk in
% this file. this is the point where the simulation starts caring about the data and
% starts calculating spectral distance values.
%
% simdata.cbinfile -  name of the cbin file for which this simdata file was
% made
%
% If use_audio_metadata =1
% It also plots a graph where each line represents a single file. the line
% plots the cdf of the difference between simulated and real distance values. 
% all files for whom the first_chunk_first_point has been calculated
% accurately, the difference value should be between -1.5 to 1.5 (although not
% strictly so). Additionally, each line is clickable. Clicking a given line
% displays the name of the file on the axes. this feature can be useful if
% you see a whacky line (with cdf vastly different than its peers), you can
% immediately get the file associated with it. 
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
%
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=5;

prob_path=pwd;
prob_path2=pwd;

in_message1='Please select the batch file';
in_message2='Please select the template file that was used in the original recording';
in_message3='Please enter the channel specification (chanspec) of the channel in the cbin file which you want to read';
in_message4='Please enter whether audio recording meta data is available for use? 0 or 1 ';
in_message5='Please enter pre-time (in seconds) that you would like to use as default in case audio metadata is not available';
in_message6='Please enter sampling frequency (in Hz) that you would like to use as default in case audio metadata is not available';

if nargin<narg_min
    [batchname,batchpath]=uigetfile([prob_path filesep 'batch*'],in_message1); 
    [chanspec]=input([in_message3 '\n-->  '],'s'); 
    [pretime]=input([in_message5 '\n-->  ']); 
    [fs]=input([in_message6 '\n-->  ']); 
    use_audio_metadata=input([in_message4 '\n-->  ']);
    if use_audio_metadata
        [templname,templpath]=uigetfile([prob_path2 filesep '*.dat'],in_message2);  
    else
       templname='';
       templpath='';
    end
end

% packaging the inputs into the inputs structure. This can be useful in
% case you need to store the inputs as meta-data with the output. 
inputs=struct('batchname',batchname,'batchpath',batchpath,'templname',templname,...
      'templpath',templpath,'chanspec',chanspec,'use_audio_metadata',...
      use_audio_metadata,'pretime',pretime);

% processing supplementary inputs

% Assigning default values to supplementary inputs
supp_inputs.window_min=50;
supp_inputs.overwrite='n';
supp_inputs.no_dists_for_compar=500;
supp_inputs.write_to_disk_q=1; % should the function write a file to disk containing its output  
supp_inputs.disk_write_dir=batchpath;


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
if ~isempty(batchpath)
    if ~strcmpi(batchpath(end),filesep)
        batchpath=[batchpath,filesep];
    end
end
if ~isempty(templpath)
    if ~strcmpi(templpath(end),filesep)
        templpath=[templpath,filesep];
    end
end

if supp_inputs.write_to_disk_q
    if ~strcmpi(supp_inputs.disk_write_dir(end),filesep)
        supp_inputs.disk_write_dir=[supp_inputs.disk_write_dir,filesep];
    end
end


% creating a figure and axes
if use_audio_metadata
    fig1=figure;
    axes1=axes('parent',fig1);
end
%

%% Body of the function

% loading the template and the batch file
fid=fopen([batchpath batchname],'r');
cbinpath=batchpath; % assumption
cbinname=fgetl(fid);

if ~use_audio_metadata
    fcfp=fix(pretime*fs); % fcfp= first chunk first point.
    if fcfp==0
        fcfp=fcfp+1;% the +1 is useful if 
    % the pre-time is set to zero. the index of fcfp then becomes 1 and not
    % zero. 
    end
else
    % loading template file
    templ=load([templpath templname]);
end

% defining prc_vec for later use
prc_vec=(0.1:0.1:100);
offset_min=0;

while ischar(cbinname)
    skip_file=0;

    [~,name,~]=fileparts(cbinname);
    simdatafile=[cbinpath name '_simdata.mat'];
    if strcmpi(supp_inputs.overwrite,'n')
        if exist(simdatafile,'file')
            warning('Not writing the *_simdata.mat file for file %s as one already exists',cbinname);
            cbinname=fgetl(fid);
            continue
        end
    end
        
    cbinname
    tic
     
    if use_audio_metadata        
        % loading rec and tmp files
        rdat=readrecf([cbinpath cbinname]); 
        if isempty(rdat)
            rec_exitsts=0;
        else
            rec_exists=1;
        end
        tmpfile=[cbinpath name '.tmp'];
        if exist(tmpfile,'file')
            real_dists=load(tmpfile);
            tmp_exists=1;
        else
            tmp_exists=0;
        end


        % extracting number of ntemplates form rdat
        no_templates=length(rdat.thresh); % this is not expected to change per file

        % verifying equality of no_templates from rdat and the template file
        if no_templates~=size(templ,2)
            h=errordlg(['The number of templates obtained from the rec file is '...
                'not equal to the number of templates from the template file. '...
                'May be the wrong template file was selected. Aborting']);
            waitfor(h)
            return
        end

        templ_len=size(templ,1);% this is not expected to change per file
        nfft=templ_len*2;% this is not expected to change per file

        % rewriting tmpdat to separate the tmvs for different templates

        tmpdat2=zeros([fix(length(real_dists)/no_templates),no_templates]);
        for ii=1:no_templates
            tmpdat2(:,ii)=real_dists(ii:no_templates:end);
        end
        real_dists=tmpdat2;
        clear tmpdat2;

        % extracting fs and pretime from rdat
        fs=rdat.adfreq;
        pretime=rdat.tbefore;

        fcfp=fix(pretime*fs); % fcfp= first chunk first point

        [dat,~]=ReadDataFile([cbinpath cbinname],chanspec);

        % initializing loop control parameters
        found_fcfp=0;
        len_buffer=10;
        new_diffs=zeros(1,2*supp_inputs.window_min+(2*len_buffer)+1);
        new_offsets=(0-supp_inputs.window_min-len_buffer:0+supp_inputs.window_min+len_buffer);
        diffs=[];
        offsets=[];

        while ~found_fcfp

            for i=1:length(new_offsets)
                relevant_dat=dat(fcfp+new_offsets(i):end);
                [sim_dists,~]=evtafsim(relevant_dat,fs,nfft,templ,0,0,0,0);

                if length(sim_dists)<supp_inputs.no_dists_for_compar
                    supp_inputs.no_dists_for_compar=length(sim_dists);
                    warning(['Using less than supp_inputs.no_dists_for_compar # of chunks for calculating the diffs for file ' cbinname '.\n The file may be too short']);
                end

                sim_dists_for_compar=sim_dists(1:supp_inputs.no_dists_for_compar,1); 
                % in case of multiple template, we will take only the distances from the first template. Any template will do for our purpose
                try
                    real_dists_for_compar=real_dists(1:supp_inputs.no_dists_for_compar,1);
                catch err
                        err
                        choice=questdlg(['Please view the command line for the error message. '... 
                                'Length of simulated tmp file and actual tmp file is not equal for file '...
                                cbinname '. One of the reasons could be a smaller than optimal values of '...
                                'supp_inputs.window_min. Skip the file or Abort the function?'],'','Skip','Abort','Skip');
                        if strcmpi(choice,'Skip')
                            skip_file=1;
                            cbinname=fgetl(fid);
                            break
                        elseif strcmpi(choice,'Abort')
                            return
                        else
                            txt=sprintf('Not writing the *_simdata.mat file for file %s',cbinname);
                            h245=warndlg(txt);
                            waitfor(h245)
                            skip_file=1;
                            cbinname=fgetl(fid);
                            break
                        end           


                end
                new_diffs(i)=sum(abs(real_dists_for_compar-sim_dists_for_compar));    
            end
            if skip_file
                break
            end

            % adding new offsets and diffs to the old ones and sorting them
            offsets=[offsets,new_offsets];
            diffs=[diffs,new_diffs];
            [offsets,indices]=sort(offsets);
            diffs=diffs(indices);

            % finding minimum difference and the its corresponding offset
            [min_diff,ind_min]=min(diffs);
            offset_min=offsets(ind_min);

            % are no. of offsets to left (or right) of minimum greater than supp_inputs.window_min
            no_offsets_to_left_enuf=ind_min-1>=supp_inputs.window_min;
            no_offsets_to_right_enuf=length(offsets)-ind_min>=supp_inputs.window_min;

            if no_offsets_to_left_enuf && no_offsets_to_right_enuf

               % the diffs to the left and right of the min_diff 
               left_diffs=diffs(ind_min-supp_inputs.window_min:ind_min-1);
               right_diffs=diffs(ind_min+1:ind_min+supp_inputs.window_min);

               if mean(left_diffs)>min_diff && mean(right_diffs)>min_diff
                   % if the min_diff lies in a valley, we have found the fcfp
                   found_fcfp=1;           
               elseif mean(left_diffs)<=min_diff && mean(right_diffs)<=min_diff
                   % expand on both sides and get new_offsets
                   new_offsets_l=(offsets(1)-supp_inputs.window_min:offsets(1)-1);
                   new_offsets_r=(offsets(end)+1:offsets(end)+supp_inputs.window_min);
                   new_offsets=[new_offsets_l,new_offsets_r];

               elseif mean(left_diffs)<min_diff
                   % expand on left side and get new_offsets
                   new_offsets=(offsets(1)-supp_inputs.window_min:offsets(1)-1);
               else 
                   % expand on the right side and get new_offsets
                   new_offsets=(offsets(end)+1:offsets(end)+supp_inputs.window_min);
               end       
            else
                if ~no_offsets_to_left_enuf
                    % expand on left and get new_offsets
                    new_offsets=(offsets(1)-supp_inputs.window_min:offsets(1)-1);
                else
                    %expand on right and get new_offsets
                    new_offsets=(offsets(end)+1:offsets(end)+supp_inputs.window_min);
                end        
            end       
        end

        if skip_file
            continue
        end
        truedat=dat(fcfp+offset_min:end);
        [sim_dists,~]=evtafsim(truedat,fs,nfft,templ,0,0,0,0);
        if size(sim_dists,1)~=size(real_dists,1)
            choice=questdlg(['Length of simulated tmp file and actual tmp file is not equal for file. '...
                cbinname 'One of the reasons could be a smaller than optimal values of supp_inputs.window_min. '...
                'Skip the file or Abort the function?'],'','Skip','Abort','Skip');

            if strcmpi(choice,'Skip')
                cbinname=fgetl(fid);
                continue
            elseif strcmpi(choice,'Abort')
                return
            else
                txt=sprintf('Not writing the *_simdata.mat file for file %s',cbinname);
                h245=warndlg(txt);
                waitfor(h245)
                cbinname=fgetl(fid);
                continue
            end           
        end

        % plotting the distribution of the differences
        true_diffs=sim_dists(:,1)-real_dists(:,1);
        percentiles=prctile(true_diffs,prc_vec);
        p=plot(axes1,percentiles,prc_vec,'ButtonDownFcn',@givefname);
        set(p,'userdata',{cbinname,axes1});
        hold(axes1,'on')
        
    end

    if supp_inputs.write_to_disk_q==1    
        simdata.first_chunk_first_point=fcfp+offset_min;
        simdata.cbinfile=cbinname;
        simdata_file=[supp_inputs.disk_write_dir name '_simdata.mat'];
        save(simdata_file,'simdata');
    end
    toc
    cbinname=fgetl(fid);
end

if exist('curr_axes','var')
  if supp_inputs.write_to_disk==1
        hold(axes1,'off')
        title(axes1,'CDF of differences between simulated and real distances for each file with first chunk first points determined by this function')
        xlabel(axes1,'Difference between distance values');
        ylabel(axes1,'Cumulative percentile')
        saveas(fig1,[supp_inputs.disk_write_dir batchname '_sim_vs_real_distances_diffs_cdf.fig'])
  end
end

% function for getting the filename after clicking on a curve on the
% figure
function givefname(gcbo,eventdata,handles)

usrdat=get(gcbo,'userdata');
fname=usrdat{1};
curr_axes=usrdat{2};
oldfname_handle=get(curr_axes,'userdata');
delete(oldfname_handle);
hte=textbp(fname,'Interpreter','None');
set(curr_axes,'userdata',hte)








