%-----------------------------------------------------------------------
% Job saved on 21-Oct-2020 11:41:40 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6685)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%% load data
clear all;
clc;
% %This part is added manually to set paths and directories before running
% %analyses
% spm_path='D:\software\study\SPM12\spm12';
% addpath(genpath(spm_path));

%directory where your files are
root_dir='D:\PhD\PPI\preprocess';
subjs={'P88','P89'};%subject filename
sessions={'block_1','block_2','block_3','block_4'};
nvolumes=167;

for i=1:length(subjs);

    subject=subjs{i};
 %% extract VOI time serious  
    spm_dir = [root_dir filesep subject filesep 'GLM\SPM.mat']; % SPM.mat in GLM file
    mask_dir = [root_dir filesep subject filesep 'GLM\mask.nii,1']; % mask.nii in GLM file
    matlabbatch{1}.spm.util.voi.spmmat = {spm_dir};
    matlabbatch{1}.spm.util.voi.adjust = NaN;
    matlabbatch{1}.spm.util.voi.session = 1;
    matlabbatch{1}.spm.util.voi.name = 'IFG'; %VOI name
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [50 30 14];%VOI peak coordinate
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3; %VOI radius
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
    matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {mask_dir}; %mask gerenated from GLM
    matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
    matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

    spm_jobman('run',matlabbatch);
    clear matlabbatch 
    
%% PPI estimate
    GLM_dir = [root_dir filesep subject filesep 'GLM'];
    matlabbatch{1}.spm.stats.ppi.spmmat = { [GLM_dir filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {[GLM_dir filesep 'VOI_IFG_1.mat']};
    matlabbatch{1}.spm.stats.ppi.type.ppi.u = [1 1 1
                                               2 1 -1];
    matlabbatch{1}.spm.stats.ppi.name = 'IFGx(Bird-Rest)';
    matlabbatch{1}.spm.stats.ppi.disp = 0;

    spm_jobman('run',matlabbatch);
    clear matlabbatch    
    
%% PPI GLM 
%% load four blocks data
    factor_onset = load([root_dir filesep subject filesep 'factor_onset.mat']);  %load onset and regressors
    load([root_dir filesep subject filesep 'GLM/PPI_IFGx(Bird-Rest).mat']); %load PPI_IFGx(Bird-Rest).mat

    files=cell(nvolumes,1); %167 volumes
    allblocks_files=[];
    for j=1:length(sessions);
        session=sessions{j};
    subj_dir= [root_dir filesep subject filesep session];
    [paths names] = filesearch_substring(subj_dir,'swafrw',0);%should give 167 files
    for k=1:nvolumes
        try
        files(k)={[paths{k} filesep names{k}]};
        catch
            disp('block missing- maybe no 3D files present?')
        end
    end
    allblocks_files=[allblocks_files;files];%should give 668 files
    end
%% Parameters1
    GLM_dir = [root_dir filesep subject filesep 'PPI'];
    matlabbatch{1}.spm.stats.fmri_spec.dir = {GLM_dir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = allblocks_files;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
%% load PPI,Block data
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'PPI-interaction';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = PPI.ppi;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'VIFG-BOLD';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = PPI.Y;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = 'Psych_Bird-Rest';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = PPI.P;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).name = 'block1';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).val = factor_onset.block1;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).name = 'block2';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).val = factor_onset.block2;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).name = 'block3';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).val = factor_onset.block3;
%% Parameters2
    motion_dir = [root_dir filesep subject filesep 'headmotion.txt'];  %% headmotion parameters path
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {motion_dir};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
%% Model estimate
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
%% Contrast
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'PPI-Interaction';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;

    spm_jobman('run',matlabbatch);
    clear matlabbatch factor_onset PPI;
    fprintf('Subject %s have done',subject);
end 
