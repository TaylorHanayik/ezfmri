function stat_1st_block (prefix, basefmriname, kTR, mocoRegress)
% first level statistics
%  basefmriname : names of 4D fMRI data (one per session)
%  kTR : repeat time in seconds
%  mocoRegress : (optional) if true include motion as regressor
%Examples
% stat_1st_blockSub('swa','fMRI.nii', 2)
% stat_1st_blockSub('swa',strvcat('fMRI1.nii','fMRI2.nii'), 2)

if ~exist('mocoRegress','var')
    mocoRegress = true; %regress motion parameters
end
kDur = 13; %duration of blocks in secs, set to 1 for events
%names for each condition
names{1} = 'leftTap';
names{2} = 'rightTap';
%c{session, cond}
onsets{1,1} = [0 52.797 105.594 158.406 211.203 264 265.094 316.797 369.594 422.406 475.203 528 529.094];
onsets{1,2} = [26.406 79.203 132 184.797 237.594 290.406 343.203 396 448.797 501.594 554.406];
%if two sessions....
%onsets{2,1} = [7.906, 39.906, 71.906, 103.906, 135.906, 167.906, 199.906, 231.906, 263.906, 295.906, 327.906];
%onsets{2,2} = [7.906, 39.906, 71.906, 103.906, 135.906, 167.906, 199.906, 231.906, 263.906, 295.906, 327.906];

% --- no need to edit below
nSessions = size(onsets,1);
nCond = size(onsets,2);
fprintf('Experiment has %d sessions with %d conditions\n',nSessions,nCond);
if (numel(names) ~= nCond)
    error('Please make sure number of condition names match number of onset rows');
end
if nSessions ~= size(basefmriname,1)
    error('There must be %d sessions of fMRI data', nSessions);
end    
%prepare SPM
if exist('spm','file')~=2; fprintf('%s requires SPM\n',which(mfilename)); return; end;
spm('Defaults','fMRI');
spm_jobman('initcfg'); % useful in SPM8 only
clear matlabbatch
%get files if not specified....
if ~exist('kTR','var') || (kTR <= 0)
    error('%s requires the repeat-time (TR) in seconds', mfilename);
end 
%next make sure each image has its full path
fmriname = [];
for ses = 1:nSessions
    [pth,nam,ext,vol] = spm_fileparts( deblank (basefmriname(ses,:)));
    if isempty(pth)
        pth = pwd; 
    end;
    fmriname = strvcat(fmriname, fullfile(pth, [nam ext])); %#ok<REMFF1>
end
%next - generate output folder
if ~exist('statdirname','var') %no input: select fMRI file[s]
 [pth,nam] = spm_fileparts(deblank(fmriname(1,:)));
 statdirname = nam;
 fprintf('Directory for SPM.mat file not specified - using folder named %s\n',statdirname);
end
predir = pwd;
%create new stat directory
if isempty(pth); pth=pwd; end;
statpth = fullfile(pth, statdirname);
if exist(statpth, 'file') ~= 7; mkdir(statpth); end;
fprintf(' SPM.mat file saved in %s\n',statdirname);
if (kDur > 5) && (kDur < 32); 
    hpf = kDur * 4;
	temporalderiv = false;
	fprintf('Block design : using %.1fs high pass filter with no temporal derivative.\n',hpf);
else
    temporalderiv = true;
    hpf = 128;
	fprintf('Event-related design : using %.1fs high pass filter with a temporal derivative.\n',hpf);
end;
% MODEL SPECIFICATION
%--------------------------------------------------------------------------
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = {statpth};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = kTR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
for ses = 1:nSessions 
    %sesFiles = fmriname{ses};%getsesvolsSubSingle(fmriname, ses);
    sesFiles = getsesvolsSubSingle(prefix, fmriname, ses);
    fprintf('Session %d has %d volumes\n',ses, length(sesFiles) );
    %sesFiles = getsesvolsSubSingle(basefmriname, ses);
    matlabbatch{1}.spm.stats.fmri_spec.sess(ses).scans = sesFiles;
    for c = 1:nCond
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).name = deblank(char(names{c}));
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).onset = cell2mat(onsets(ses, c));
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).duration = kDur;
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
    end;
    matlabbatch{1}.spm.stats.fmri_spec.sess(ses).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(ses).regress = struct('name', {}, 'val', {});
    if mocoRegress
        [p,n] = spm_fileparts(deblank(fmriname(ses,:)));
        motionFile = fullfile(p, ['rp_', n, '.txt']);
        if ~exist(motionFile, 'file')
            error('Unable to find realign parameters for motion correction %s', motionFile);
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).multi_reg = {motionFile};
    else
        matlabbatch{1}.spm.stats.fmri_spec.sess(ses).multi_reg = {''}; 
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess(ses).hpf = hpf;
end
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
if temporalderiv
	matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
else 
	matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
end;
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
% MODEL ESTIMATION
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(statpth,'SPM.mat'));
% INFERENCE
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.con.spmmat = cellstr(fullfile(statpth,'SPM.mat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name    = 'Task>Rest';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = ones(1,nCond);
nContrast = 1;
if nCond > 1
    for pos = 1: nCond
        for neg = 1: nCond
            if pos ~= neg
                nContrast = nContrast + 1;
                c = zeros(1,nCond);
                c(pos) = 1;
                c(neg) = -1;
                matlabbatch{3}.spm.stats.con.consess{nContrast}.tcon.convec = c;
                matlabbatch{3}.spm.stats.con.consess{nContrast}.tcon.name = [char(names{pos}) '>' char(names{neg})];
            end % j ~= i
        end %for j
    end %for i
end % > 1 conditions
if temporalderiv %zero pad temporal derivations
    for c = 1 : nContrast
        for cond = 1 : nCond
            v = matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec;
            c2 = cond * 2;
            v = [v(1:(c2-1)) 0 v(c2:end)];
            matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec = v;
            
        end
    end
end
if mocoRegress %add 6 nuisance regressors for motion paramets (rotation + translation]
    for c = 1 : nContrast
    	 matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec = [ matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec  0 0 0 0 0 0];
    end  
end
if (nSessions > 1) %replicate contrasts for each session
    for c = 1 : nContrast
        v = matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec;
        for s = 2 : nSessions
         matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec = [matlabbatch{3}.spm.stats.con.consess{c}.tcon.convec v];
        end
    end
end
spm_jobman('run',matlabbatch);
cd(predir); %return to starting directory...
%stat_1st_blockSub

%%%%  LOCAL SUB FUNCTIONS FOLLOW

function [sesvols] = getsesvolsSubSingle(prefix, fmriname, session)
%* Load all fMRI images from single sessions
[pth,nam,ext,vol] = spm_fileparts( deblank (fmriname(session,:))); %#ok<*NASGU>
sesname = fullfile(pth,[prefix, nam, ext]);
hdr = spm_vol(sesname);
nvol = length(hdr);
if (nvol < 2), fprintf('Error 4D fMRI data required %s\n', sesname); return; end;
sesvols = cellstr([sesname,',1']);
for vol = 2 : nvol
    sesvols = [sesvols; [sesname,',',int2str(vol)]  ]; %#ok<AGROW>
end;
%end getsesvolsSubSingle()