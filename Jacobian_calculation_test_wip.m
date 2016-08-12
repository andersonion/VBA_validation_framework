%% Data Simulator for VBM Analysis
%  Purpose: to manually induce known variations in a given image or set of
%  images.
%
%addpath('/home/rja20/cluster_code/workstation_code/shared/mathworks/NIfTI_20140122/');
AntsPath = '/cm/shared/apps/ANTS/';

%%  Define Variables

%dir_name = '/glusterspace/BJ/jacobian_calculation_testing_2/';
dir_name = '/glusterspace/BJ/JacTest_for_paper/'

%runnos = {'mcnamara'};
runno = 'mcnamara';
%runno = 'obrien';

%runnos = {'PEDS_t1'};

%test_erosion = 1;

concentric = 1;
%ants_geo_option = 0;
fix_DCM = 0;
%use_inverse_warp = 1;
morph_with_ants = 0;

v_field_smooth_radius = 1; % 0 is the fastest, but it shows up in the Jacobian. On the bright side, it might be easier to detect bad calcs this way.
% Change the above back to 1 for full on data analysis

d_radii = (-6:1:6);
d_radii(d_radii == 0)=[];

test_erosion=(d_radii < 1);


for nn = 1:length(d_radii)%n = 1:length(runnos)
    
    if test_erosion(nn)
        direction_string = 'eroded';
        ants_morph_command = 'ME';
    else
        direction_string = 'dilated';
        ants_morph_command = 'MD';
    end
    
    % if use_inverse_warp
    %     warp_string = 'inverse';
    % else
    %     warp_string = 'forward';
    % end
    
    %runno = runnos{n}
    test_nii_path =[dir_name runno '.nii.gz']; % functionize
    
    radius = 30; % Radius of original (moving) test sphere
    %d_radius = 3; % Change in radius of fixed test sphere
    d_radius = d_radii(nn)
    
    reg_mask_tolerance = 4; % A larger mask will be created for registration
    
    if test_erosion(nn)
        new_radius = radius + d_radius;
        big_radius = radius + reg_mask_tolerance;
    else % test_dilation
        new_radius = radius + d_radius;
        big_radius = new_radius + reg_mask_tolerance;
    end
    
    id_string = ['_r' num2str(radius) '_' direction_string '_to_r' num2str(new_radius)];
    
    orig_mask_path = [dir_name runno '_radius' num2str(radius)  '_mask.nii.gz'];
    new_mask_path = [dir_name runno '_radius' num2str(new_radius) '_mask.nii.gz'];
    big_mask_path = [dir_name runno '_radius' num2str(big_radius) '_mask.nii.gz'];
    
    if ~concentric
        geo_string = 'tangential';
    else
        geo_string = 'concentric';
    end
    
    warp_path = [dir_name runno id_string '_' ];
    
    %% Load data
    master_nii=load_untouch_nii(test_nii_path);
    
    if fix_DCM
        master_nii.hdr.hist.quatern_d=1;%%%
        master_nii.hdr.hist.srow_x(1)=-1*master_nii.hdr.hist.srow_x(1);%%%
        master_nii.hdr.hist.srow_y(2)=-1*master_nii.hdr.hist.srow_y(2);%%%
        %master_nii.hdr.hist.srow_z(3)=-1*master_nii.hdr.hist.srow_z(3);% Unusual option...testing to find a case where ANTs works
        master_nii.hdr.hist.srow_z(4)=-1*master_nii.hdr.hist.srow_z(4);%%%
        master_nii.hdr.hist.qoffset_z(1)=-1*master_nii.hdr.hist.qoffset_z(1);%%%
    end
    
    %% Make morphological changes
    sizes = size(master_nii.img);
    largest_dim = find(sizes == max(sizes));
    largest_dim = largest_dim(1); % In case of two dims equally being the largest
    
    shift_val = 20;%round(min(sizes)/5);%10; % Offset the sphere to ensure "robustness to assymetry"
    
    extra = zeros(size(sizes));
    extra(largest_dim) = d_radius;
    
    start=-1*floor(sizes/2);
    ending = start+sizes-1;
    gv_x = (start(1):1:ending(1)) + shift_val;
    gv_y = (start(2):1:ending(2)) + shift_val;
    gv_z = (start(3):1:ending(3)) + shift_val;
    
    [X,Y,Z] = meshgrid(gv_y,gv_x,gv_z);
    
    if ~concentric
        [X2,Y2,Z2] = meshgrid((gv_y+extra(2)),(gv_x+extra(1)),(gv_z+extra(3)));
    else
        X2=X;
        Y2=Y;
        Z2=Z;
    end
    
    
    r2 = radius*radius;
    sphere_array = ((X.*X + Y.*Y + Z.*Z)<=r2);
    
    if ~morph_with_ants
        new_r2 = new_radius*new_radius;
        new_sphere_array = ((X2.*X2 + Y2.*Y2 + Z2.*Z2)<=new_r2);
    end
    
    big_r2 = big_radius*big_radius;
    big_sphere_array = ((X.*X + Y.*Y + Z.*Z)<=big_r2);
    
    
    if ~exist(orig_mask_path,'file')
        orig_mask = master_nii;
        orig_mask.img = sphere_array;
        save_untouch_nii(orig_mask,orig_mask_path);
    end
    
    if ~exist(new_mask_path,'file')
        if morph_with_ants
            erode_cmd = [AntsPath 'ImageMath 3 ' new_mask_path ' ' ants_morph_command ' ' orig_mask_path ' ' num2str(abs(d_radius))]
            system(erode_cmd)
            %if ~exist(big_mask_path,'file')
            %    new_mask = load_untouch_nii(new_d_path);
            %    new_sphere = new_mask.img
            %end
        else
            new_mask = master_nii;
            new_mask.img = new_sphere_array;
            save_untouch_nii(new_mask,new_mask_path)
        end
    end
    
    if ~exist(big_mask_path,'file')
        big_mask = master_nii;
        big_mask.img = big_sphere_array;
        save_untouch_nii(big_mask,big_mask_path);
    end
    
    fixed_reg_params = '-c [3000x3000x3000x3000,1.e-7,20] -s 4x2x1x0.5vox -f 6x4x2x1 -l 1 -u 1 -z 1'; % No full sampling for prototyping--turn for production!
    
    forward_warp = [warp_path '0Warp.nii.gz'];
    
    for use_inverse_warp = [1,0];
        
        if use_inverse_warp
            warp_string = 'inverse'
            warp_file = [warp_path '0InverseWarp.nii.gz'];
            inverse_index = 2;
        else
            warp_string = 'forward';
            warp_file = [warp_path '0Warp.nii.gz'];
            inverse_index = 1;
        end
        
        if ~exist(warp_file,'file')
            reg_cmd = [AntsPath 'antsRegistration -d 3 -m MeanSquares[' new_mask_path ',' orig_mask_path ',1,4] -x [' big_mask_path ',' big_mask_path '] -t SyN[0.4,1,' num2str(v_field_smooth_radius) '] ' fixed_reg_params ' -o ' warp_path]
            system(reg_cmd)
        end
        
        
        %% Apply the warps to the original masks for directly calculating volume change
        final_mask_path = [dir_name runno id_string '_warped.nii.gz'];
        if ~exist(final_mask_path,'file')
            apply_cmd = [AntsPath 'antsApplyTransforms --float -d 3 -i ' orig_mask_path ' -o ' final_mask_path ' -t ' forward_warp ' -r ' orig_mask_path ' -n NearestNeighbor']
            system(apply_cmd)
        end
        
        % Calculate forward or inverse Jacobian images
        geo_log_jacobian_path = [dir_name runno id_string '_' warp_string '_geo_logJacobian.nii.gz'];
        if ~exist(geo_log_jacobian_path,'file')
            geo_log_jacobian_cmd = [AntsPath 'CreateJacobianDeterminantImage 3 ' warp_file ' ' geo_log_jacobian_path ' 1 1'] % num2str(ants_geo_option)
            system(geo_log_jacobian_cmd)
        end
        
        FD_log_jacobian_path = [dir_name runno id_string '_' warp_string '_FD_logJacobian.nii.gz'];
        if ~exist(FD_log_jacobian_path,'file')
            FD_log_jacobian_cmd = [AntsPath 'CreateJacobianDeterminantImage 3 ' warp_file ' ' FD_log_jacobian_path ' 1 0'] % num2str(ants_geo_option)
            system(FD_log_jacobian_cmd)
        end
        
        geo_jacobian_path = [dir_name runno id_string '_' warp_string '_geo_Jacobian.nii.gz'];
        if ~exist(geo_jacobian_path,'file')
            geo_jacobian_cmd = [AntsPath 'CreateJacobianDeterminantImage 3 ' warp_file ' ' geo_jacobian_path ' 0 1']; %  num2str(ants_geo_option)
            system(geo_jacobian_cmd)
        end
        
        FD_jacobian_path = [dir_name runno id_string '_' warp_string '_FD_Jacobian.nii.gz'];
        if ~exist(FD_jacobian_path,'file')
            FD_jacobian_cmd = [AntsPath 'CreateJacobianDeterminantImage 3 ' warp_file ' ' FD_jacobian_path ' 0 0' ];% num2str(ants_geo_option)
            system(FD_jacobian_cmd)
        end
        
        
        
        geo_jac_nii=load_untouch_nii(geo_jacobian_path);
        geo_jac_data = geo_jac_nii.img;
        
        FD_jac_nii=load_untouch_nii(FD_jacobian_path);
        FD_jac_data = FD_jac_nii.img;
        
        geo_log_jac_nii=load_untouch_nii(geo_log_jacobian_path);
        geo_log_jac_data = geo_log_jac_nii.img;
        
        FD_log_jac_nii=load_untouch_nii(FD_log_jacobian_path);
        FD_log_jac_data = FD_log_jac_nii.img;
        
        
        orig_mask_nii = load_untouch_nii(orig_mask_path);
        orig_mask = orig_mask_nii.img;
        
        new_mask_nii = load_untouch_nii(new_mask_path);
        new_mask = new_mask_nii.img;
        
        warped_mask_nii = load_untouch_nii(final_mask_path);
        warped_mask = warped_mask_nii.img;
        warped_mask_vox_volume = sum(warped_mask(:));
        
        orig_indices = find(orig_mask);
        original_vox_volume = length(orig_indices);
        
        new_indices = find(new_mask);
        new_vox_volume = length(new_indices);
        
        new_geo_jacs = geo_jac_data(orig_indices);
        new_geo_vox_volume_1e = sum(new_geo_jacs(:));
        measured_geo_Jacobian(nn,inverse_index) = new_geo_vox_volume_1e/original_vox_volume
        
        new_FD_jacs = FD_jac_data(orig_indices);
        new_FD_vox_volume_1e = sum(new_FD_jacs(:));
        measured_FD_Jacobian(nn,inverse_index) = new_FD_vox_volume_1e/original_vox_volume
        
        
        %mean_geo_log_jac(n) = mean(geo_log_jac_data(orig_indices));
        average_geo_log_jac(nn,inverse_index) = mean(geo_log_jac_data(orig_indices));
        
        converted_geo_log_jac(nn,inverse_index) = exp(average_geo_log_jac(nn,inverse_index));
        
        average_FD_log_jac(nn,inverse_index) = mean(FD_log_jac_data(orig_indices));
        
        converted_FD_log_jac(nn,inverse_index) = exp(average_FD_log_jac(nn,inverse_index));
    end
    
    %expected_change_in_volume(n) = new_vox_volume/original_vox_volume
    actual_Jacobian(nn,inverse_index) = warped_mask_vox_volume/original_vox_volume
end
