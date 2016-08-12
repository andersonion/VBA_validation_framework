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

plot_stuff = 0;

master_dir = '/glusterspace/BJ/Mcnamara_labels';
runnos = {'64944' '64953' '64959' '64962' '64968' '64974' '65394' '65408' '65411' '65414'}

jac_array = [];


for rr = 1:length(runnos)
    runno =['S' runnos{rr}];
    work_dir = '/glusterspace/BJ/phantom_inputs_eval/';
    
    optional_suffix = ['_' runno];
    
    jac_labels = [23 24];
    labels = jac_labels;
    
    label_names = {'CPu' 'Hc'};% 'ac' 'fi' 'cc' 'ic'};
    local_labels = [master_dir runno '_labels_in_Control_space.nii.gz'];
    
    %%
    r_dilate = 2;
    if ~already_analyzed && ~shells_exist
        label_nii=load_untouch_nii(local_labels);
        label_data = label_nii.img;
        
        group_mask = ismember(label_data,jac_labels);
        
        for l_index = 1:length(labels)
            
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
        for cc = 1
            contrast = contrasts{cc};
            effects_image = [master_dir runno '_f_log_jacobian.nii'];
            if ~exist(effects_image,'file')
                effects_image = [effects_image '.gz'];
            end
            
            effects_nii=load_untouch_nii(effects_image);
            effects = effects_nii.img;
            
            relevant_labels = jac_labels;
            % effects = 10.^(effects); %unlogging
            
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
            
            jac_array = cat(1,jac_array,array);
            
        end
    end
    
end
file_3 = [work_dir 'phantom_jac_results'];
if already_analyzed
    load(file_3)
else
    save(file_3,'jac_array');
end

for cc = 1
    
    current_data = jac_array;
    num_labels = 2;
    if plot_stuff
        figure(117*cc + num_labels)
    end
    
    
    for LL = 1:num_labels
        
        mean_effect = current_data(LL:num_labels:end,2);
        std_effect = current_data(LL:num_labels:end,3);
        mean_leak = current_data(LL:num_labels:end,4);
        std_leak = current_data(LL:num_labels:end,5);
        %snr = mean_effect./std_leak;
        %snr = mean_effect./mean_leak;
        %snr = (mean_effect.*std_leak)./(mean_leak.*std_effect);
        snr = (mean_effect-mean_leak)./((0.5*((std_effect).^2+(std_leak).^2)).^(1/2));% Misnomer--actually the sensitivity index
        if plot_stuff
            for pp = 1:3
                index = ((LL-1)*3+pp);
                h1 = subplot(num_labels,3,index);
                if pp == 1
                    data_to_plot = mean_effect;
                    bar(data_to_plot,'LineWidth',1.5)
                    hold on
                    h3 = errorbar((1:10),data_to_plot,std_effect,'k','LineStyle','none','LineWidth',1.5);
                elseif pp ==2
                    data_to_plot = mean_leak;
                    bar(data_to_plot,'LineWidth',1.5)
                    hold on
                    h3 = errorbar((1:10),data_to_plot,std_leak,'k','LineStyle','none','LineWidth',1.5);
                else
                    data_to_plot = snr;
                    bar(data_to_plot,'LineWidth',1.5)
                    hold on
                end
                
                colors = cat(1,hsv(10),[0 0 0],hsv(10),[0 0 0]);
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
end