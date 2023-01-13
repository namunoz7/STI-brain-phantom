function CreateSimulatedData(ModelParams,SeqParams,SimParams)


% ModelParams.R1map_file
% ModelParams.R2starmap_file
% ModelParams.M0map_file
% ModelParams.Segmentation_file
% ModelParams.Chimap_file
% ModelParams.Segmentation_file
% ModelParams.BrainMask_file

% SeqParams.TR          Repetition time in secs
% SeqParams.TE          Echo time in secs
% SeqParams.FlipAngle   flip angle in degrees

% SimParams.PhaseOffset multiplier term of a quadratic phase over the brain
%                       0 no phase offset; 2pi phase difference inside
%                       brain mask
% SimParams.Shimm       boolean 0 if no additional shimms are applied 1 if
%                       2nd order shimms are applied inside brain mask
% SimParams.BrainExtracted  boolean 0 makes the simulation with whole head
%                       susceptibility, 1 makes only head
% SimParams.Res         Resolution of output
% SimParams.Output_dir  Output Directory

% modeln refers to which susceptibility model is being used
% BackGroundFieldRemovalToSimulate refers to whether there is or not
% background field removal
% SimulateData is wheather the simulation part is bing done or not



%% Data needed to be loaded for simulation
% clear
% SimulateData =1 ;
% modeln = [1 2 3 4]
% BackGroundFieldRemovalToSimulate=[1 3]
% %    BackGroundFieldRemoval = 3 ; % no background field removal
% %    BackGroundFieldRemoval = 1 ; % removes only the outside  the brain region
% %    BackGroundFieldRemoval = 2 ; % removes also the calcification


Sufix(1).name='BrainExtracted';
Sufix(2).name='';
savedata = 1;

% some relevant parameters for a simulation
B0=7;
B0_dir=[0 0 1];
gyro = 42.57747892;% MHz/T



% baseDir = '/project/3015069.02/HeadPhantomData/7T/Derived/sandbox/'
% cd(baseDir)
% a=genpath([baseDir,'func/']);
% addpath(a)
% a=genpath('/home/mrphys/josmar/MATLAB/MathWorkDownloads/reisert-unring-ddd43da65219')
% addpath(a)
modelname = SimParams.Output_dir ;

Folder=[modelname,'/GroundTruthHR/'];

%%
    for BackGroundFieldRemoval = [1 2]

        
        
        M0=load_nii(ModelParams.M0map_file);
        R1=load_nii(ModelParams.R1map_file);
        R2star=load_nii(ModelParams.R2starmap_file);
        chi0=load_nii(ModelParams.Chimap_file);
        FinalSegment=load_nii(ModelParams.Segmentation_file);
        Brain=load_nii(ModelParams.BrainMask_file);
        Brain=Brain.img;
        voxel_size = round(M0.hdr.dime.pixdim(2:4)*100)/100;
        
        
        chi = chi0 ;
        
        
        % 'BrainExtracted';
        if BackGroundFieldRemoval == 1
            chi.img = (chi.img - mean(chi.img(Brain == 1))) .* Brain ;
            %             fv(Brain)
            %             fv(chi.img)
            Shimm = 0 ; % there is never a Bo shimming being applied when the the simulation is brain only
        end
        
        if BackGroundFieldRemoval == 2
            Shimm = SimParams.Shimm ;
        end
        %% creates dipole Kernel
        dims = size(M0.img);
        %         D = create_dipole_kernel([0 0 1], voxel_size, dims, 1);
        %         field= real(ifftn(fftn(chi.img).*(D)))  ;% if chi is in ppm than field is in ppm
        %       doing some padding
        D = create_dipole_kernel(B0_dir , voxel_size, 2 * dims, 1);
        
        chitemp = ones( 2 * dims) * chi.img(end,end,end);
        chitemp (1:dims(1),1:dims(2),1:dims(3)) = chi.img;
        field= real(ifftn(fftn(chitemp).*(D)))  ;% if chi is in ppm than field is in ppm
        
        field = field(1:dims(1),1:dims(2),1:dims(3));
        clear chitemp
        clear D
        
        
        %% Brain shimming
        
        if Shimm==0
            fieldb.img=field;
        else
            [~,fieldb,~]=PolyFitShimLike(make_nii(field),make_nii(single(Brain)),2);
        end;
        %% Brain Phase Offset
        
        if SimParams.PhaseOffset==0 % 
            PhaseOffset = 0;
        else
            [c , w ] = centerofmass(M0.img)
            [y,x,z] = meshgrid([1:dims(2)]-c(2), [1:dims(1)]-c(1), [1:dims(3)]-c(3));
            temp = (x/w(1)).^2 + (y/w(2)).^2 + (z/w(3)).^2 ;
            PhaseOffset = - temp/(max(temp(Brain==1))-min(temp(Brain==1)))*pi*SimParams.PhaseOffset;
        end
        
        %%
        
        try
            figureJ(1),set(gcf,'Position',[1,29,1920,977])
            subplot(221)
            Orthoview(Brain)
            %    title('shim region')
            %    subplot(222)
            %   Orthoview(fieldb.img)
            subplot(222)
            Orthoview(fieldb.img,[],[-1 1]*2.5)
            title('shimmed field')
            %    subplot(224)
            %    Orthoview(field)
            subplot(224)
            Orthoview(field,[],[-1 1]*2.5)
            title('unshimmed field')
            
            subplot(223)
            
            
            Orthoview(angle(exp(i*PhaseOffset)));
            title('B1 Phase Offset')
            if savedata == 1
                savefig([Folder,'ShimmingAndB1phase',Sufix(BackGroundFieldRemoval).name])
            end
        catch
            display('error displaying results')
        end;
        field = fieldb.img;
        clear fieldb
        
        
        %%
        %             for res = SeqParams.res
        
        for TE = SeqParams.TE
            
            if savedata == 1
                %                 if model==1
                FolderData=[modelname,'/Simulated_',num2str(floor(mean(SimParams.Res))),'p',num2str(floor(10*(mean(SimParams.Res)-floor(mean(SimParams.Res))))),'mm/']
                FolderHR=[modelname,'/SimulatedHR/']
                %                 else
                %                     FolderData=['ChallengeDataOverModulated/Simulated_',num2str(floor(res)),'p',num2str(floor(10*(res-floor(res)))),'mm/']
                %                     FolderHR=['ChallengeDataOverModulated/SimulatedHR/']
                %
                %                 end;
                
                if exist(FolderData,'dir')~=7
                    mkdir(FolderData)
                end
                
                if exist(FolderHR,'dir')~=7
                    mkdir(FolderHR)
                end
                
            end
            % config is a structure with all sort of parameters
            SequenceParam.TR=SeqParams.TR;                     %repetition time (in seconds);
            SequenceParam.TE=TE;                     %echo time (in seconds);
            SequenceParam.theta=SeqParams.FlipAngle;                     %flip angle in degrees;
            if length(SimParams.Res)==1
                SequenceParam.res=[1 1 1 ]*SimParams.Res;             %resolution in mm;
            elseif length(SimParams.Res)==3
                SequenceParam.res=SimParams.Res;             %resolution in mm;
            else
                display('you have not entered the resolution in the correct format, we [ 1 1 1 ] will be assumed')
                SequenceParam.res=[1 1 1 ];             %resolution in mm;
            end
            TissueParam.M0=double(M0.img);              %Water Concentration (arbitrry Units);
            TissueParam.R1=double(R1.img);                      %longitudinal relaxation rate (1/s);
            TissueParam.R2star=R2star.img;                  %apparent transverse relaxation rate (1/s);
            TissueParam.field=field * B0 * gyro;                   %field should be (1/s); - this assumes the field was calculated in ppm
            %     this to have the info of the
            TissueParam.PhaseOffset = PhaseOffset;             %PhaseOfsett at TE=0;
            TissueParam.res = voxel_size       ;             %resolution in mm;
            
            [sig,vox_out,sigHR]= DataSimulation(SequenceParam,TissueParam);
            %% checking the complex data created at the two resolutions
            slicepercent=40 ; % this is only to create a data visualization
            dims_start=size(sigHR);
            dims_end=size(sig);
            figureJ(3)
            set(gcf,'Position',[1,29,1920,977])
            colormap(gray)
            subplot(2,2,1)
            imab(abs(sig(:,:,round(slicepercent*dims_end(3)/100)))),title('Low Res Magnitude')
            subplot(2,2,2)
            imab(angle(sig(:,:,round(slicepercent*dims_end(3)/100)))),title('Low Res Phase')
            
            subplot(2,2,3)
            imab(abs(sigHR(:,:,round(slicepercent*dims_start(3)/100)))),title('High Res Magnitude')
            subplot(2,2,4)
            imab(angle(sigHR(:,:,round(slicepercent*dims_start(3)/100)))),title('High Res Phase')
            
            if savedata == 1
                savefig([FolderData,'HighvsLowRes_TE',num2str(TE*1000),Sufix(BackGroundFieldRemoval).name])
                save_nii(make_nii(abs(sig),vox_out),[FolderData,'Magnitude_TE',num2str(TE*1000,'%03.f'),Sufix(BackGroundFieldRemoval).name,'.nii.gz'])
                save_nii(make_nii(angle(sig),vox_out),[FolderData,'Phase_TE',num2str(TE*1000,'%03.f'),Sufix(BackGroundFieldRemoval).name,'.nii.gz'])
            end
            if savedata == 1
                save_nii(make_nii(abs(sigHR),voxel_size),[FolderHR,'Magnitude_TE',num2str(TE*1000,'%03.f'),Sufix(BackGroundFieldRemoval).name,'.nii.gz'])
                save_nii(make_nii(angle(sigHR),voxel_size),[FolderHR,'Phase_TE',num2str(TE*1000,'%03.f'),Sufix(BackGroundFieldRemoval).name,'.nii.gz'])
            end
            figureJ(3)
            set(gcf,'Position',[1,29,1920,977])
            colormap(gray)
            subplot(2,2,1)
            Orthoview(abs(sig)),title('Low Res Magnitude')
            subplot(2,2,2)
            Orthoview(angle(sig)),title('Low Res Phase')
            
            subplot(2,2,3)
            Orthoview(abs(sigHR)),title('High Res Magnitude')
            subplot(2,2,4)
            Orthoview(angle(sigHR)),title('High Res Phase')
            
        end
        
        
        %%
        %% checking the lowres susceptibility ground truth
        dims_start=size(sigHR);
        dims_end=size(sig);
        clear sig
        figureJ(5)
        set(gcf,'Position',[1,29,1920,977])
        % probably more correct - but does not match the data signal
        % simulation part
        % chi_end = real( ifftshift(ifftn(ifftshift(crop(fftshift(fftn(fftshift(chi.img.*Brain))),dims_end)))))*prod(dims_end)/prod(dims_start);
        % probably less correct - but does match the data signal
        % simulation part
         chi_end = real( (ifftn(ifftshift(crop(fftshift(fftn((chi.img))),dims_end)))))*prod(dims_end)/prod(dims_start);
        chi_end = permute(unring(permute(unring(chi_end),[3 1 2])),[2 3 1])	;
        [X,Y,Z] = ndgrid(single(1:dims_start(1)),single(1:dims_start(2)),single(1:dims_start(3)));
        [Xq,Yq,Zq] = ndgrid(single(linspace(1,dims_start(1),dims_end(1))),...
            single(linspace(1,dims_start(2),dims_end(2))),single(linspace(1,dims_start(3),dims_end(3))));
        
        chi_end_interp = interpn(X,Y,Z,single(chi.img.*Brain),Xq,Yq,Zq);
        
        
        
        subplot(311)
        Orthoview(chi.img.*Brain,[],[-0.04 0.08]);colorbar; % high res susceptibility map
        title('High Res susceptibility')
        subplot(312)
        Orthoview(chi_end,[],[-0.04 0.08]);colorbar; % high res susceptibility map
        title('Low Res susceptibility FFTcrop')
        subplot(313)
        Orthoview(chi_end_interp,[],[-0.04 0.08]);colorbar; % high res susceptibility map
        title('Low Res susceptibility Interp')
        if savedata == 1
            savefig([FolderData,'GroundTruthSusceptibility',Sufix(BackGroundFieldRemoval).name])
            %             save_nii(make_nii(Brain,vox_out),[FolderData,'Brain',Sufix(BackGroundFieldRemoval).name,'.nii.gz'])
            save_nii(make_nii(chi_end_interp,vox_out),[FolderData,'Chi_interp',Sufix(BackGroundFieldRemoval).name,'.nii.gz'])
            save_nii(make_nii(chi_end,vox_out),[FolderData,'Chi_crop',Sufix(BackGroundFieldRemoval).name,'.nii.gz'])
            %                 end
            clear chi_end
            clear chi_end_interp
            %% testing the downsampling of the segmentation model
            % I would have probably prefered to do this segmentation  on the
            % probability maps used to create the ground truth
            % susceptibility...
            temp=zeros([size(Xq), max(FinalSegment.img(:))]);
            for label=1:max(FinalSegment.img(:))
                
                temp(:,:,:,label) = interpn(X,Y,Z,single(FinalSegment.img==label),Xq,Yq,Zq);
                %                 makes the computation really slow
                %                 tmp = real( (ifftn(ifftshift(crop(fftshift(fftn(single(FinalSegment.img==label))),dims_end)))))*prod(dims_end)/prod(dims_start);
                %                 temp(:,:,:,label) = permute(unring(permute(unring(tmp),[3 1 2])),[2 3 1])	;
                
                
            end
            [temp_val temp_pos]= max(temp,[],4);
            Segment_LR=temp_pos;
            subplot(211),Orthoview(Segment_LR)
            subplot(212),Orthoview(FinalSegment.img)
            temp= interpn(X,Y,Z,single(Brain),Xq,Yq,Zq);
            Brain_LR= temp>0.9;
            if savedata == 1
                save_nii(make_nii(Segment_LR,vox_out),[FolderData,'FinalSegment',Sufix(BackGroundFieldRemoval).name,'.nii.gz'])
                save_nii(make_nii(single(Brain_LR),vox_out),[FolderData,'Brain',Sufix(BackGroundFieldRemoval).name,'.nii.gz'])
            end
            clear temp, clear temp_val, clear temp_pos, clear Segment_LR, clear Brain_LR
            
            %%
            
            % Checking how it compares to real data that is stored in
            % data/raw
            %
            try
                TE = [3 11.5 20 28.5]*10e-3;
                DataDir='data/raw/';
                files{1}.name='MP2RAGEME_INV2';
                files{2}.name='MP2RAGEME_INV2_e2';
                files{3}.name='MP2RAGEME_INV2_e3';
                files{4}.name='MP2RAGEME_INV2_e4';
                [temp, pos]= min(SequenceParam.TE-TE);
                MagnData = load_nii([DataDir,files{pos}.name,'.nii']);
                PhaseData = load_nii([DataDir,files{pos}.name,'ph.nii']);
                figureJ(10)
                set(gcf,'Position',[1,29,1920,977])
                subplot(2,2,1)
                imab((MagnData.img(:,:,round(slicepercent*dims_start(3)/100)))),title('Real Magnitude')
                subplot(2,2,2)
                imab(PhaseData.img(:,:,round(slicepercent*dims_start(3)/100))),title('Real Phase')
                subplot(2,2,3)
                imab(abs(sigHR(:,:,round(slicepercent*dims_start(3)/100)))),title('High Res Magnitude')
                subplot(2,2,4)
                imab(angle(sigHR(:,:,round(slicepercent*dims_start(3)/100)))),title('High Res Phase')
                
                subplot(2,2,1)
                Orthoview((MagnData.img)),title('Real Magnitude')
                subplot(2,2,2)
                Orthoview(PhaseData.img),title('Real Phase')
                subplot(2,2,3)
                Orthoview(abs(sigHR)),title('High Res Magnitude')
                subplot(2,2,4)
                Orthoview(angle(sigHR)),title('High Res Phase')
                clear MagnData
                clear PhaseData
                clear sigHR
                if savedata == 1
                    savefig([Folder,'ComparisonToRealData',Sufix(BackGroundFieldRemoval).name])
                end
            end;
        end;
        
    end;
    save( [ modelname, '/SimulationParameters.mat'],'ModelParams','SeqParams','SimParams')
 