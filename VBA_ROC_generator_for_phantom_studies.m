%% A script for analyzing phantom VBA results

%   Define:  Label(s) of interest
%           MDT Label set
%           Contrasts of interest
%           Effect size image (from VBA results - should be impervious to
%               direction)
%


%   Load MDT Labels
%   Write results

%%
AntsPath = '/cm/shared/apps/ANTS/';

already_analyzed = 0;
shells_exist = 1;


master_dir = '/glusterspace/VBM_13mcnamara02_DTI101b_SyN_1_3_1-work/';
%SyNs = {'p5_3_1' '1_1_0' '1_3_1' '5_3_1' '5_3_3'};
SyNs = {'p5_1_0' 'p5_3_1' 'p5_3_3' '1_1_0' '1_3_1' '1_3_3' '5_1_0' '5_3_1' '5_3_3'};

if 0
    comparison_string = 'Control_gt_Phantom';
    alt_comp_string = 'Control_gt_vs';
    comparison_string_2 = 'Phantom_gt_Control';
    alt_comp_string_2 = 'Phantom_gt_vs';
else
    comparison_string = 'Phantom_gt_Control';
    alt_comp_string = 'Phantom_gt_vs';
    comparison_string_2 = 'Control_gt_Phantom';
    alt_comp_string_2 = 'Control_gt_vs';
end


jac_array = [];

for SS = 1:length(SyNs)
    SyN = SyNs{SS};
    %work_dir = '/glusterspace/BJ/ROCs_for_phantom_surfstat/';
    work_dir = '/glusterspace/BJ/reanalyzed_ROCs_for_phantom_surfstat/';
    SyN_string =  ['SyN_' SyN];
    optional_suffix = ['_' SyN_string];
    
    jac_labels = [23 24];
    labels =  [jac_labels] ;
    label_names = {'CPu' 'Hc'};
    
    % Should we use the 1_3_1 labelset, as that is what defined each
    % subject's ROIs?
    %MDT_labelset = [master_dir 'dwi/' SyN_string '_fa/faMDT_Control_n10a/median_images/MDT_labels_DTI101b.nii.gz' ];
    MDT_labelset = [master_dir 'dwi/' SyN_string '_fa/faMDT_Control_n10d/median_images/MDT_labels_DTI101b.nii.gz' ];
    MDT_mask = [master_dir 'dwi/' SyN_string '_fa/faMDT_Control_n10d/median_images/MDT_mask_e3.nii' ];
    contrasts = {'jac'};
    %vba_results_path = [master_dir 'dwi/' SyN_string '_fa/faMDT_Control_n10a/vbm_analysis/faMDT_Control_n10a_3vox_smoothing-results/surfstat/'];
    vba_results_path = [master_dir 'dwi/' SyN_string '_fa/faMDT_Control_n10d/vbm_analysis/faMDT_Control_n10d_3vox_smoothing-results/surfstat/'];
    %vba_results_path = [master_dir 'dwi/' SyN_string '_fa/faMDT_Control_n10d/vbm_analysis/faMDT_Control_n10d_3vox_smoothing-results/antsr/'];
    
    
    %%
    
    if ~already_analyzed
        label_nii=load_untouch_nii(MDT_labelset);
        label_data = label_nii.img;
        CPu_mask = ismember(label_data,jac_labels(1));
        Hc_mask = ismember(label_data,jac_labels(2));
        
        mask_nii= load_untouch_nii(MDT_mask);
        full_mask = mask_nii.img;
        
        size_mask = length(find(full_mask));
        size_true_CPu = length(find(CPu_mask));
        size_true_Hc = length(find(Hc_mask));
        size_false_CPu = size_mask - size_true_CPu;
        size_false_Hc = size_mask - size_true_Hc;
        
        
        for cc = 1%length(contrasts)
            contrast = contrasts{cc};
            
            p_image_CPu = [vba_results_path contrast '/' contrast '_pvalunc_' comparison_string '.nii'];
            
            
            if ~exist(p_image_CPu,'file')
                p_image_CPu = [p_image_CPu '.gz'];
            end
            
            if ~exist(p_image_CPu,'file')
                p_image_CPu = [vba_results_path contrast '/' contrast '_pvalunc_' alt_comp_string '.nii'];
                
            end
            
            if ~exist(p_image_CPu,'file')
                p_image_CPu = [p_image_CPu '.gz'];
            end
            
            p_nii_CPu=load_untouch_nii(p_image_CPu);
            p_vals_CPu = p_nii_CPu.img;
            
            
            q_image_CPu = [vba_results_path contrast '/' contrast '_qval_' comparison_string '.nii'];
            
            
            if ~exist(q_image_CPu,'file')
                q_image_CPu = [q_image_CPu '.gz'];
            end
            
            if ~exist(q_image_CPu,'file')
                q_image_CPu = [vba_results_path contrast '/' contrast '_qval_' alt_comp_string '.nii'];
                %effects_image = [vba_results_path contrast '/' contrast '_effect_' alt_comp_string '.nii'];
            end
            
            if ~exist(q_image_CPu,'file')
                q_image_CPu = [q_image_CPu '.gz'];
            end
            
            q_nii_CPu=load_untouch_nii(q_image_CPu);
            q_vals_CPu = q_nii_CPu.img;
            
            
            p_image_Hc = [vba_results_path contrast '/' contrast '_pvalunc_' comparison_string_2 '.nii'];
            
            
            if ~exist(p_image_Hc,'file')
                p_image_Hc = [p_image_Hc '.gz'];
            end
            
            if ~exist(p_image_Hc,'file')
                p_image_Hc = [vba_results_path contrast '/' contrast '_pvalunc_' alt_comp_string_2 '.nii'];
            end
            
            if ~exist(p_image_Hc,'file')
                p_image_Hc = [p_image_Hc '.gz'];
            end
            
            p_nii_Hc=load_untouch_nii(p_image_Hc);
            p_vals_Hc = p_nii_Hc.img;
            
            
            q_image_Hc = [vba_results_path contrast '/' contrast '_qval_' comparison_string_2 '.nii'];
            
            
            if ~exist(q_image_Hc,'file')
                q_image_Hc = [q_image_Hc '.gz'];
            end
            
            if ~exist(q_image_Hc,'file')
                q_image_Hc = [vba_results_path contrast '/' contrast '_qval_' alt_comp_string_2 '.nii'];
                %effects_image = [vba_results_path contrast '/' contrast '_effect_' alt_comp_string '.nii'];
            end
            
            if ~exist(q_image_Hc,'file')
                q_image_Hc = [q_image_Hc '.gz'];
            end
            
            q_nii_Hc=load_untouch_nii(q_image_Hc);
            q_vals_Hc = q_nii_Hc.img;
            
            
            
            false_region_CPu = find(full_mask.*(1-CPu_mask));
            false_region_Hc = find(full_mask.*(1-Hc_mask));
            true_region_CPu = find(CPu_mask);
            true_region_Hc = find(Hc_mask);
            
            false_pvals_CPu = p_vals_CPu(false_region_CPu);
            true_pvals_CPu = p_vals_CPu(true_region_CPu);
            
            false_qvals_CPu = q_vals_CPu(false_region_CPu);
            true_qvals_CPu = q_vals_CPu(true_region_CPu);
            
            false_pvals_Hc = p_vals_Hc(false_region_Hc);
            true_pvals_Hc = p_vals_Hc(true_region_Hc);
            
            false_qvals_Hc = q_vals_Hc(false_region_Hc);
            true_qvals_Hc = q_vals_Hc(true_region_Hc);
            
            
            th_step = 0.005;
            th_step_2 = 0.0005;
            t_point = 0.99;
            index = 1;
            for th=[(0:th_step:t_point) ((t_point+th_step_2):th_step_2:1)]
                
                
                true_and_false_p_CPu(SS,index,2) = length(find(false_pvals_CPu <= th))/size_false_CPu ;
                true_and_false_p_CPu(SS,index,1) = length(find(true_pvals_CPu <= th))/size_true_CPu;
                
                true_and_false_q_CPu(SS,index,2) = length(find(false_qvals_CPu <= th))/size_false_CPu ;
                true_and_false_q_CPu(SS,index,1) = length(find(true_qvals_CPu <= th))/size_true_CPu;
                
                %
                %
              
                true_and_false_p_Hc(SS,index,2) = length(find(false_pvals_Hc <= th))/size_false_Hc ;
                true_and_false_p_Hc(SS,index,1) = length(find(true_pvals_Hc <= th))/size_true_Hc;
                
                true_and_false_q_Hc(SS,index,2) = length(find(false_qvals_Hc <= th))/size_false_Hc ;
                true_and_false_q_Hc(SS,index,1) = length(find(true_qvals_Hc <= th))/size_true_Hc;
                
                index = index+1;
            end
            
            
        end
    end
    
end
file_1 = [work_dir 'ROC_data_p_CPu'];
file_2 = [work_dir 'ROC_data_q_CPu'];
file_3 = [work_dir 'ROC_data_p_Hc'];
file_4 = [work_dir 'ROC_data_q_Hc'];
if already_analyzed
    load(file_1)
    load(file_2)
    load(file_3)
    load(file_4)
else
    save(file_1,'true_and_false_p_CPu');
    save(file_2,'true_and_false_q_CPu');
    save(file_3,'true_and_false_p_Hc');
    save(file_4,'true_and_false_q_Hc');
end

colors = hsv(9);

hizzle=figure(11)
set(hizzle,'Color',[1 1 1])
for ss = 1:SS
    plot(true_and_false_p_CPu(ss,:,2),true_and_false_p_CPu(ss,:,1),'Color',colors(ss,:),'LineWidth',2)
    hold on
    set(gca,'FontName','Ariel','FontSize',14,'FontWeight','bold');
    set(gca,'XTick', [0,0.2,0.4,0.6,0.8,1],'XTickLabel', [0,0.2,0.4,0.6,0.8,1]);
    set(gca,'YTick', [0,0.2,0.4,0.6,0.8,1],'YTickLabel', [0,0.2,0.4,0.6,0.8,1]);
end

hizzle=figure(12)
set(hizzle,'Color',[1 1 1])
for ss = 1:SS
    plot(true_and_false_q_CPu(ss,:,2),true_and_false_q_CPu(ss,:,1),'Color',colors(ss,:),'LineWidth',2)
    hold on
    set(gca,'FontName','Ariel','FontSize',14,'FontWeight','bold');
    set(gca,'XTick', [0,0.2,0.4,0.6,0.8,1],'XTickLabel', [0,0.2,0.4,0.6,0.8,1]);
    set(gca,'YTick', [0,0.2,0.4,0.6,0.8,1],'YTickLabel', [0,0.2,0.4,0.6,0.8,1]);
end

hizzle=figure(21)
set(hizzle,'Color',[1 1 1])
for ss = 1:SS
    plot(true_and_false_p_Hc(ss,:,2),true_and_false_p_Hc(ss,:,1),'Color',colors(ss,:),'LineWidth',2)
    hold on
    set(gca,'FontName','Ariel','FontSize',14,'FontWeight','bold');
    set(gca,'XTick', [0,0.2,0.4,0.6,0.8,1],'XTickLabel', [0,0.2,0.4,0.6,0.8,1]);
    set(gca,'YTick', [0,0.2,0.4,0.6,0.8,1],'YTickLabel', [0,0.2,0.4,0.6,0.8,1]);
end

hizzle=figure(22)
set(hizzle,'Color',[1 1 1])
for ss = 1:SS
    plot(true_and_false_q_Hc(ss,:,2),true_and_false_q_Hc(ss,:,1),'Color',colors(ss,:),'LineWidth',2)
    hold on
    set(gca,'FontName','Ariel','FontSize',14,'FontWeight','bold');
    set(gca,'XTick', [0,0.2,0.4,0.6,0.8,1],'XTickLabel', [0,0.2,0.4,0.6,0.8,1]);
    set(gca,'YTick', [0,0.2,0.4,0.6,0.8,1],'YTickLabel', [0,0.2,0.4,0.6,0.8,1]);
end



if 0
    for LL = 1:num_labels
        
        mean_effect = current_data(LL:num_labels:end,2);
        std_effect = current_data(LL:num_labels:end,3);
        mean_leak = current_data(LL:num_labels:end,4);
        std_leak = current_data(LL:num_labels:end,5);
        %snr = mean_effect./std_leak;
        %snr = mean_effect./mean_leak;
        %snr = (mean_effect.*std_leak)./(mean_leak.*std_effect);
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
