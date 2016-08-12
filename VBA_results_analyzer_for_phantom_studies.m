%% A script for analyzing phantom VBA results

%   Define:  Label(s) of interest
%           MDT Label set
%           Contrasts of interest
%           Effect size image (from VBA results - should be impervious to
%               direction)
%


%   Load MDT Labels
%       create label subset
%       dilate label subset (2 voxels)
%       subtract any and all original label regions from dilated regions
%           these dilated shells will be used for assessing "leakage"
%   Load Effect size image for contrast(n)
%       For each Label
%       calculate average effect size in ROI
%       calculate stdev of effect size in ROI
%
%       calculate average effect size (leakage)
%       calculate stdev of effect size (leakage)
%   Write results

%%
AntsPath = '/cm/shared/apps/ANTS/';

already_analyzed = 1;
shells_exist = 1;

surfstat = 1 ; % else ANTsR
log_jac = 0;
master_dir = '/glusterspace/VBM_13mcnamara02_DTI101b_SyN_1_3_1-work/';

SyNs = {'p5_1_0' 'p5_3_1' 'p5_3_3' '1_1_0' '1_3_1' '1_3_3' '5_1_0' '5_3_1' '5_3_3'};
%comparison_string = 'Control_gt_Phantoms';
%alt_comp_string = 'Control_gt_vs';
comparison_string = 'Phantom_gt_Control'; % I believe this contrast will make the sign of the effect size consistent with atrophy/expansion
alt_comp_string = 'Phantom_gt_vs'; % In an earlier iteration it was "Phantoms" not "Phantom"

jac_array = [];
fa_array = [];
dwi_array = [];

for SS = 1:length(SyNs)
    SyN = SyNs{SS};
    %work_dir = '/glusterspace/BJ/phantom_and_optimization_eval_good_ANTsR/';
    work_dir = '/glusterspace/BJ/Mcnamara_reanalyzed/';
    SyN_string =  ['SyN_' SyN];
    optional_suffix = ['_' SyN_string];
    
    jac_labels = [23 24];
    fa_labels = [10 25 26 16];
    labels =  [jac_labels fa_labels] ;
    label_names = {'CPu' 'Hc' 'ac' 'fi' 'cc' 'ic'};
    
    % Should we use the 1_3_1 labelset, as that is what defined each
    % subject's ROIs?
    %MDT_labelset = [master_dir 'dwi/' SyN_string '_fa/faMDT_Control_n10a/median_images/MDT_labels_DTI101b.nii.gz' ];
    MDT_labelset = [master_dir 'dwi/' SyN_string '_fa/faMDT_Control_n10d/median_images/MDT_labels_DTI101b.nii.gz' ];
    contrasts = {'jac' 'fa' 'dwi'};
    %vba_results_path = [master_dir 'dwi/' SyN_string '_fa/faMDT_Control_n10a/vbm_analysis/faMDT_Control_n10a_3vox_smoothing-results/surfstat/'];
    if surfstat
        vba_results_path = [master_dir 'dwi/' SyN_string '_fa/faMDT_Control_n10d/vbm_analysis/faMDT_Control_n10d_3vox_smoothing-results/surfstat/'];
    else
        vba_results_path = [master_dir 'dwi/' SyN_string '_fa/faMDT_Control_n10d/vbm_analysis/faMDT_Control_n10d_3vox_smoothing-results/antsr/'];
    end
    
    %%
    r_dilate = 2;
    if ~already_analyzed && ~shells_exist
        label_nii=load_untouch_nii(MDT_labelset);
        label_data = label_nii.img;
        
        fa_mask = ismember(label_data,fa_labels);
        jac_mask = ismember(label_data,jac_labels);
        
        for l_index = 1:length(labels)
            
            if l_index < 3
                group_mask = jac_mask;
            else
                group_mask = fa_mask;
            end
            label = labels(l_index);
            label_mask = ismember(label_data,label);
            label_mask_path = [work_dir 'label_' num2str(label)  optional_suffix '_mask.nii.gz'];
            
            if ~exist(label_mask_path,'file')
                label_mask_nii = label_nii;
                label_mask_nii.img = label_mask;
                save_untouch_nii(label_mask_nii,label_mask_path);
            end
            
            dilated_label_path = [work_dir 'r' num2str(r_dilate) '_dilated_mask_for_label_' num2str(label) optional_suffix '.nii.gz'];
            
            if ~exist(dilated_label_path,'file')
                dilate_cmd = [AntsPath 'ImageMath 3 ' dilated_label_path ' MD ' label_mask_path ' ' num2str(r_dilate)]
                system(dilate_cmd)
            end
            
            dilated_nii=load_untouch_nii(dilated_label_path);
            dilated_data = dilated_nii.img;
            
            leakage_shell = dilated_data - group_mask;
            leakage_shell(leakage_shell < 0) = 0;
            
            leakage_mask_path = [work_dir 'r' num2str(r_dilate) '_leakage_mask_for_label_' num2str(label) optional_suffix '.nii.gz'];
            
            if ~exist(leakage_mask_path,'file')
                leakage_mask_nii = label_nii;
                leakage_mask_nii.img = leakage_shell;
                save_untouch_nii(leakage_mask_nii,leakage_mask_path);
            end
        end
    end
    
    if ~already_analyzed
        for cc = 1%length(contrasts)
            contrast = contrasts{cc};
            if surfstat
                effects_image = [vba_results_path contrast '/' contrast '_effect_' comparison_string '.nii'];
            else
                effects_image = [vba_results_path contrast '/' comparison_string '_effect.nii'];
            end
            
            if ~exist(effects_image,'file')
                effects_image = [effects_image '.gz'];
            end
            
            if ~exist(effects_image,'file')
                effects_image = [vba_results_path contrast '/' contrast '_effect_' alt_comp_string '.nii'];
            end
            
            if ~exist(effects_image,'file')
                effects_image = [effects_image '.gz'];
            end
            effects_image
            effects_nii=load_untouch_nii(effects_image);
            effects = effects_nii.img;
            
            if strcmp(contrast,'jac')
                relevant_labels = jac_labels;
                if ~log_jac
                %%%%effects = 10.^(effects); %unlogging --Le Whoops!
                    effects = exp(effects);
                end
            elseif strcmp(contrast,'fa')
                relevant_labels = fa_labels;
                %relevant_labels=[];
            else
                relevant_labels = labels;
                %relevant_labels=[] ;    
            end
            
            array = zeros([length(relevant_labels) 5]);
            
            for l_index = 1:length(relevant_labels)
                label = relevant_labels(l_index);
                label_mask_path = [work_dir 'label_' num2str(label)  optional_suffix '_mask.nii.gz'];
                label_nii=load_untouch_nii(label_mask_path);
                label_mask = label_nii.img;
                
                leakage_mask_path = [work_dir 'r' num2str(r_dilate) '_leakage_mask_for_label_' num2str(label) optional_suffix '.nii.gz'];
                leakage_nii=load_untouch_nii(leakage_mask_path);
                leakage_mask = leakage_nii.img;
                
                effect_vector = effects(label_mask==1);
                leakage_vector = effects(leakage_mask==1);
                
                array(l_index,1) = label;
                array(l_index,2) = mean(effect_vector(:));
                array(l_index,3) = std(effect_vector(:));
                array(l_index,4) = mean(leakage_vector(:));
                array(l_index,5) = std(leakage_vector(:));
            end
            
            
            if strcmp(contrast,'jac')
                jac_array = cat(1,jac_array,array);
            elseif strcmp(contrast,'fa')
                fa_array = cat(1,fa_array,array);
            else
                dwi_array = cat(1,dwi_array,array);
            end
            
        end
    end

end
    file_1 = [work_dir 'phantom_fa_results'];
    file_2 = [work_dir 'phantom_dwi_results'];
    if log_jac
         file_3 = [work_dir 'phantom_jac_results'];
    else
        file_3 = [work_dir 'phantom_unlogged_jac_results'];
    end
if already_analyzed
    load(file_1)
    load(file_2)
    load(file_3)
else
 %   save(file_1,'fa_array');
 %   save(file_2,'dwi_array');
    save(file_3,'jac_array');
end

for cc = 1:3
    if cc ==1
        current_data = dwi_array;
        num_labels = 6;
    elseif cc ==2
        current_data = fa_array;
        num_labels = 4;
    elseif cc==3
        current_data = jac_array;
        num_labels = 2;
    end
    figure(10*cc + num_labels)
    
   
     
    for LL = 1:num_labels

        mean_effect = current_data(LL:num_labels:end,2);
        std_effect = current_data(LL:num_labels:end,3);
        mean_leak = current_data(LL:num_labels:end,4);
        std_leak = current_data(LL:num_labels:end,5);
        snr = (mean_effect-mean_leak)./((0.5*((std_effect).^2+(std_leak).^2)).^(1/2));% Misnomer--actually the sensitivity index
    
        for pp = 1:3
            index = ((LL-1)*3+pp);
            h1 = subplot(num_labels,3,index);
            if pp == 1
                data_to_plot = mean_effect;
                bar(data_to_plot,'LineWidth',1.5)
                hold on
                h3 = errorbar((1:9),data_to_plot,std_effect,'k','LineStyle','none','LineWidth',1.5);
            elseif pp ==2
                data_to_plot = mean_leak;
                bar(data_to_plot,'LineWidth',1.5)
                hold on
                h3 = errorbar((1:9),data_to_plot,std_leak,'k','LineStyle','none','LineWidth',1.5);
            else
                data_to_plot = snr;
                bar(data_to_plot,'LineWidth',1.5)
                hold on
            end
            
            colors = cat(1,hsv(9),[0 0 0],hsv(9),[0 0 0]);
            for i = 1:numel(data_to_plot)
                bar(i, data_to_plot(i), 'parent', h1, 'facecolor', colors(i,:));
            end
            if pp ==1
                %[lower_lim,upper_lim] = get(h1,'YLim');
                limits = get(h1,'YLim');
                set(h1,'YLim',limits);
            elseif pp==2
                set(h1,'YLim',limits);
            end
            
            set(h1,'XTickLabel','','XTick', []);
            set(h1,'TickLength',[0,0]);
            %set(h1, 'pos', pp(index,:));
            set(h1, 'LineWidth', 1.5);
            set(h1,'FontName','Ariel','FontSize',16,'FontWeight','bold');
            hold off
        end
    end
end