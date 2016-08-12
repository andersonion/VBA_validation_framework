%% Data Simulator for VBM Analysis
%  Purpose: to manually induce known variations in a given image or set of
%  images.  Current changes include FA and cortical dialation.
%
%addpath('/home/rja20/cluster_code/workstation_code/shared/mathworks/NIfTI_20140122/');
AntsPath = '/cm/shared/apps/ANTS/';

%%  Define Variables

dir_name = '/glusterspace/BJ/jacobian_calculation_testing_2/';

runnos = {'hismom'};

matlab_jacs = 1;
rotate_around_z = 1;
fix_DCM = 1;

for n = 1:length(runnos)
    runno = runnos{n}

    test_nii_path =[runno '.nii.gz']; % functionize
    
    orig_e_path = [dir_name runno '_orig_eroded_mask.nii.gz'];
    new_e_path = [dir_name runno '_eroded_mask.nii.gz'];
    
    orig_d_path = [dir_name runno '_orig_dilate_mask.nii.gz'];
    new_d_path = [dir_name runno '_dilated_mask.nii.gz'];
    
    orig_mask_path = [dir_name runno '_original_mask.nii.gz'];
    new_mask_path = [dir_name runno '_new_mask.nii.gz'];
    
    warp_path = [dir_name runno '_shapeshift_'];
    
    %%
    
    r_erode = 1;
    r_dilate = 1; % Size, in voxels of dilation radius.
    
    erode_labels = [24]; % Functional label values = (NOMINAL - 1)
    
    
    dilate_labels = [23]; % Functional label values = (NOMINAL - 1)
    
    %% Load data
    
    
    orig_labels = [dir_name test_nii_path];
    master_nii=load_untouch_nii(orig_labels);
    %if fix_DCM
        master_nii.hdr.hist.quatern_d=1;%%%
        master_nii.hdr.hist.srow_x(1)=-1*master_nii.hdr.hist.srow_x(1);%%%
        master_nii.hdr.hist.srow_y(2)=-1*master_nii.hdr.hist.srow_y(2);%%%
        master_nii.hdr.hist.srow_z(4)=-1*master_nii.hdr.hist.srow_z(4);%%%
        master_nii.hdr.hist.qoffset_z(1)=-1*master_nii.hdr.hist.qoffset_z(1);%%%
    %end
    
    label_data = master_nii.img;
    
    
    %% Make morphological changes
    sizes = size(label_data);
    size_x = sizes(1)
    size_y = sizes(2)
    size_z = sizes(3)
    min_size =min(sizes);
    
    base_r = 30;
    %dil_r = base_r + r_dilate;
    %ero_r = base_r - r_erode;
    shft = 10;
    os = min_size/2-shft; %offset
    os2 = min_size/2 + shft;
    label_data(:) = 0;
    if copulate_with_cortex
        for xx=1:size_x
            for yy = 1:size_y
                for zz = 1:size_z
                    if (((xx-os)*(xx-os)+(yy-os)*(yy-os)+(zz-os)*(zz-os)) <= base_r*base_r);
                        label_data(xx,yy,zz) = erode_labels(1);
                    end
                    
                    if ((xx-os2)*(xx-os2)+(yy-os2)*(yy-os2)+(zz-os2)*(zz-os2)) <= base_r*base_r;
                        %    label_data(xx,yy,zz) = dilate_labels(1);
                    end
                    
                end
            end
        end
        
        
        
        if ~exist(orig_e_path,'file')
            erode_mask = ismember(label_data,erode_labels);
            e_mask = master_nii;
            e_mask.img = erode_mask;
            if fix_DCM
                e_mask.hdr.hist.quatern_d=1;%%%
                e_mask.hdr.hist.srow_x(1)=-1*e_mask.hdr.hist.srow_x(1);%%%
                e_mask.hdr.hist.srow_y(2)=-1*e_mask.hdr.hist.srow_y(2);%%%
                e_mask.hdr.hist.srow_z(4)=-1*e_mask.hdr.hist.srow_z(4);%%%
                e_mask.hdr.hist.qoffset_z(1)=-1*e_mask.hdr.hist.qoffset_z(1);%%%
            end
            save_untouch_nii(e_mask,orig_e_path)
        end
        
        if ~exist(orig_d_path,'file')
            dilate_mask = ismember(label_data,dilate_labels);
            d_mask = master_nii;
            d_mask.img = dilate_mask;
            if fix_DCM
                d_mask.hdr.hist.quatern_d=1;%%%
                d_mask.hdr.hist.srow_x(1)=-1*d_mask.hdr.hist.srow_x(1);%%%
                d_mask.hdr.hist.srow_y(2)=-1*d_mask.hdr.hist.srow_y(2);%%%
                d_mask.hdr.hist.srow_z(4)=-1*d_mask.hdr.hist.srow_z(4);%%%
                d_mask.hdr.hist.qoffset_z(1)=-1*d_mask.hdr.hist.qoffset_z(1);%%%
            end
            save_untouch_nii(d_mask,orig_d_path);
        end
        
        if ~exist(orig_mask_path,'file')
            orig_mask = master_nii;
            orig_mask.img = (dilate_mask + erode_mask);
            if fix_DCM
                orig_mask.hdr.hist.quatern_d=1;%%%
                orig_mask.hdr.hist.srow_x(1)=-1*orig_mask.hdr.hist.srow_x(1);%%%
                orig_mask.hdr.hist.srow_y(2)=-1*orig_mask.hdr.hist.srow_y(2);%%%
                orig_mask.hdr.hist.srow_z(4)=-1*orig_mask.hdr.hist.srow_z(4);%%%
                orig_mask.hdr.hist.qoffset_z(1)=-1*orig_mask.hdr.hist.qoffset_z(1);%%%
            end
            save_untouch_nii(orig_mask,orig_mask_path);
        end
        
        
        
        if ~exist(new_e_path,'file')
            erode_cmd = [AntsPath 'ImageMath 3 ' new_e_path ' ME ' orig_e_path ' ' num2str(r_erode)]
            system(erode_cmd)
        end
        
        if ~exist(new_d_path,'file')
            dilate_cmd = [AntsPath 'ImageMath 3 ' new_d_path ' MD ' orig_d_path ' ' num2str(r_dilate)]
            system(dilate_cmd)
        end
        
        dilated_nii=load_untouch_nii(new_d_path);
        eroded_nii=load_untouch_nii(new_e_path);
        
        dilated_data = dilated_nii.img;
        eroded_data = eroded_nii.img;
        
        big_mask_path = [dir_name runno '_big_mask.nii'];
        if ~exist(big_mask_path,'file')
            big_mask = master_nii;
            big_mask.img = ((dilated_data + erode_mask) > 0);
            if fix_DCM
                big_mask.hdr.hist.quatern_d=1;%%%
                big_mask.hdr.hist.srow_x(1)=-1*big_mask.hdr.hist.srow_x(1);%%%
                big_mask.hdr.hist.srow_y(2)=-1*big_mask.hdr.hist.srow_y(2);%%%
                big_mask.hdr.hist.srow_z(4)=-1*big_mask.hdr.hist.srow_z(4);%%%
                big_mask.hdr.hist.qoffset_z(1)=-1*big_mask.hdr.hist.qoffset_z(1);%%%
            end
            save_untouch_nii(big_mask,big_mask_path);
            
            
            dilate_mask_cmd = [AntsPath 'ImageMath 3 ' big_mask_path ' MD ' big_mask_path ' 4']
            system(dilate_mask_cmd)
        end
        
        if ~exist(new_mask_path,'file')
            new_mask_data = ((dilated_data + eroded_data) > 0);
            new_mask = master_nii;
            new_mask.img = new_mask_data;
            new_mask.hdr.hist.quatern_d=1;%%%
            new_mask.hdr.hist.srow_x(1)=-1*new_mask.hdr.hist.srow_x(1);%%%
            new_mask.hdr.hist.srow_y(2)=-1*new_mask.hdr.hist.srow_y(2);%%%
            new_mask.hdr.hist.srow_z(4)=-1*new_mask.hdr.hist.srow_z(4);%%%
            new_mask.hdr.hist.qoffset_z(1)=-1*new_mask.hdr.hist.qoffset_z(1);%%%
            save_untouch_nii(new_mask,new_mask_path);
        end
        
        
        if ~exist([warp_path '1Warp.nii.gz'],'file')
            %reg_cmd = [AntsPath 'antsRegistration -d 3 -r [' new_mask_path ',' orig_mask_path ',1] -m MeanSquares[' new_mask_path ',' orig_mask_path ',1,4] -x [' big_mask_path ',' big_mask_path '] -c [3000x3000x3000x0,1.e-8,20] -t SyN[1,3,3] -s 4x2x1x0.5vox -f 6x4x2x1 -l 1 -u 1 -z 1 -o ' warp_path]
            reg_cmd = [AntsPath 'antsRegistration -d 3 -r [' new_mask_path ',' orig_mask_path ',1] -m MeanSquares[' new_mask_path ',' orig_mask_path ',1,4] -x [' big_mask_path ',' big_mask_path '] -c [3000x3000x3000x0,1.e-8,20] -t SyN[0.4,1,0] -s 4x2x1x0.5vox -f 6x4x2x1 -l 1 -u 1 -z 1 -o ' warp_path]
            system(reg_cmd)  % Original used 'CC' metric out of habit...trying 'MeanSquares' based on http://nipy.org/dipy/examples_built/syn_registration_2d.html
        end
        
        
        %final_mask_path = [dir_name runno '_warped_mask.nii'];
        %if ~exist(final_mask_path,'file')
        %    apply_cmd = [AntsPath 'antsApplyTransforms --float -d 3 -i ' orig_mask_path ' -o ' final_mask_path ' -t ' warp_path '1Warp.nii.gz -r ' orig_mask_path ' -n NearestNeighbor'] % originally -n NearestNeighbor --but Why?
        %    system(apply_cmd)
        %end
        
        final_d_mask_path = [dir_name runno '_warped_dilated_mask.nii.gz'];
        if ~exist(final_d_mask_path,'file')
            apply_cmd = [AntsPath 'antsApplyTransforms --float -d 3 -i ' orig_d_path ' -o ' final_d_mask_path ' -t ' warp_path '1Warp.nii.gz -r ' orig_d_path ' -n NearestNeighbor'] % originally -n NearestNeighbor --but Why?
            system(apply_cmd)
        end
        
        final_e_mask_path = [dir_name runno '_warped_eroded_mask.nii.gz'];
        if ~exist(final_e_mask_path,'file')
            apply_cmd = [AntsPath 'antsApplyTransforms --float -d 3 -i ' orig_e_path ' -o ' final_e_mask_path ' -t ' warp_path '1Warp.nii.gz -r ' orig_e_path ' -n NearestNeighbor'] % originally -n NearestNeighbor --but Why?
            system(apply_cmd)
        end
        
        f_jacobian_path = [dir_name runno '_f_log_jacobian.nii'];
        if ~exist(f_jacobian_path,'file')
            in_warp = [warp_path '1Warp.nii.gz'];
            if matlab_jacs
                calculate_jacobian(in_warp,f_jacobian_path,1,rotate_around_z);
            else
                f_jacobian_cmd = [AntsPath 'CreateJacobianDeterminantImage 3 ' in_warp ' ' f_jacobian_path ' 1 0']
                system(f_jacobian_cmd)
            end
        end
        % f_no_geo_jacobian_path = [dir_name runno '_f_log_no_geo_jacobian.nii'];
        % if ~exist(f_no_geo_jacobian_path,'file')
        %     f_no_geo_jacobian_cmd = [AntsPath 'CreateJacobianDeterminantImage 3 ' warp_path '1Warp.nii.gz ' f_no_geo_jacobian_path ' 1 0']
        %     system(f_no_geo_jacobian_cmd)
        % end
        
        f_actual_jacobian_path = [dir_name runno '_f_actual_jacobian.nii'];
        if ~exist(f_actual_jacobian_path,'file')
            in_warp = [warp_path '1Warp.nii.gz'];
            if matlab_jacs
                calculate_jacobian(in_warp,f_actual_jacobian_path,0,rotate_around_z);
            else
                f_actual_jacobian_cmd = [AntsPath 'CreateJacobianDeterminantImage 3 ' in_warp ' ' f_actual_jacobian_path ' 0 1']
                %f_actual_jacobian_cmd = [AntsPath 'CreateJacobianDeterminantImage 3 ' in_warp ' ' f_actual_jacobian_path ' 0 0']
                system(f_actual_jacobian_cmd)
            end
        end
        
        % f_no_geo_actual_jacobian_path = [dir_name runno '_f_no_geo_jacobian.nii'];
        % if ~exist(f_no_geo_actual_jacobian_path,'file')
        %     f_no_geo_actual_jacobian_cmd = [AntsPath 'CreateJacobianDeterminantImage 3 ' warp_path '1Warp.nii.gz ' f_no_geo_actual_jacobian_path ' 0 0']
        %     system(f_no_geo_actual_jacobian_cmd)
        % end
        
        
        
        orig_dilated_mask_nii = load_untouch_nii(orig_d_path);
        orig_dilated_mask = orig_dilated_mask_nii.img;
        
        orig_eroded_mask_nii = load_untouch_nii(orig_e_path);
        orig_eroded_mask = orig_eroded_mask_nii.img;
        
        
        new_dilated_mask_nii = load_untouch_nii(new_d_path);
        new_dilated_mask = new_dilated_mask_nii.img;
        
        new_eroded_mask_nii = load_untouch_nii(new_e_path);
        new_eroded_mask = new_eroded_mask_nii.img;
        
        
        f_jac_nii=load_untouch_nii(f_actual_jacobian_path);
        f_jac_data = f_jac_nii.img;
        
        f_log_jac_nii=load_untouch_nii(f_jacobian_path);
        f_log_jac_data = f_log_jac_nii.img;
        
        %         if ~matlab_jacs
        %             f_no_geo_jac_nii=load_untouch_nii(f_no_geo_actual_jacobian_path);
        %             f_no_geo_jac_data = f_no_geo_jac_nii.img;
        %         end
        %
        warped_d_mask_nii = load_untouch_nii(final_d_mask_path);
        warped_d_mask = warped_d_mask_nii.img;
        warped_d_mask_vox_volume = sum(warped_d_mask(:));
        
        warped_e_mask_nii = load_untouch_nii(final_e_mask_path);
        warped_e_mask = warped_e_mask_nii.img;
        warped_e_mask_vox_volume = sum(warped_e_mask(:));
        
        
        orig_dilated_indices = find(orig_dilated_mask);
        original_d_vox_volume = length(orig_dilated_indices);
        
        orig_eroded_indices = find(orig_eroded_mask);
        original_e_vox_volume = length(orig_eroded_indices);
        
        new_dilated_indices = find(new_dilated_mask);
        new_d_vox_volume = length(new_dilated_indices);
        
        new_eroded_indices = find(new_eroded_mask);
        new_e_vox_volume = length(new_eroded_indices);
        
        dilated_jacs = f_jac_data(orig_dilated_indices);
        new_vox_volume_1d = sum(dilated_jacs(:));
        geo_measured_dilated_change_in_volume(n) = new_vox_volume_1d/original_d_vox_volume
        mean_d_log_jac(n) = mean(f_log_jac_data(orig_dilated_indices))
        u_mean_d_log_jac(n) = 10^mean_d_log_jac(n)
        
        eroded_jacs = f_jac_data(orig_eroded_indices);
        new_vox_volume_1e = sum(eroded_jacs(:));
        geo_measured_eroded_change_in_volume(n) = new_vox_volume_1e/original_e_vox_volume
        mean_e_log_jac(n) = mean(f_log_jac_data(orig_eroded_indices))
        u_mean_e_log_jac(n) = 10^mean_e_log_jac(n)
        
        %         if ~matlab_jacs
        %             dilated_no_geo_jacs = f_no_geo_jac_data(orig_dilated_indices);
        %             new_vox_volume_2d = sum(dilated_no_geo_jacs(:));
        %             no_geo_measured_dilated_change_in_volume(n) = new_vox_volume_2d/original_d_vox_volume
        %
        %             eroded_no_geo_jacs = f_no_geo_jac_data(orig_eroded_indices);
        %             new_vox_volume_2e = sum(eroded_no_geo_jacs(:));
        %             no_geo_measured_eroded_change_in_volume(n) = new_vox_volume_2e/original_e_vox_volume
        %
        %         end
        
        expected_eroded_change_in_volume(n) = new_e_vox_volume/original_e_vox_volume
        expected_dilated_change_in_volume(n) = new_d_vox_volume/original_d_vox_volume
        
        actual_dilated_change_in_volume(n) = warped_d_mask_vox_volume/original_d_vox_volume
        actual_eroded_change_in_volume(n) = warped_e_mask_vox_volume/original_e_vox_volume
    end
end
