
% The following code needs pulseq path
fov = [224e-3 224e-3 192e-3];     % Define FOV
res = [1,1,1] * 1e-3;             % voxel size

N = round(fov ./ res);

Nx = N(1); 
Ny = N(2); 
Nz = N(3);            % Define FOV and resolution

bandwidth = 200;
flip_angle = 15; % degrees

use_wave_y = 0;
use_wave_z = 0;

Tread = 1 / bandwidth;
Tpre= 3.0e-3;
riseTime = 300e-6;
Ndummy = 100;
poidisc = 1;
blip = 1;

%% Poisson-disc sampling pattern generation
bart_path='../bart';
addpath([bart_path '/matlab']);
cur_dir=pwd;
cd(bart_path);
setenv('TOOLBOX_PATH', pwd);
cd(cur_dir);

mask = ones(Ny,Nz);


if (poidisc == 0)
    R = 1;
end

% R12

mask = squeeze(readcfl('poidisc_mask'));


if (blip == 1)
    params.blip_shift = (squeeze(readcfl('table2')));
end

%% Some Hardware Constraints
 %prisma
B0field=2.89;% 6.98
sys = mr.opts('MaxGrad',60,'GradUnit','mT/m',...
    'MaxSlew',160,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0',B0field);  %max slew rate: 200

seq=mr.Sequence(sys);           % Create a new sequence object

%% design sequence 
%--------------------------------------------------------------------------

params.NcyclesGz = 8;                  % num sine cycles
params.gZmax = 1.8/1.5;                     % in Gauss/cm (i.e. 18 mT/m)
params.sZmax = 16000;                   % in Gauss/cm/s (i.e. 160 T/m/s)

params.NcyclesGy = 8;                  % num cosine cycles
params.gYmax = 1.8/1.5;                     % in Gauss/cm
params.sYmax = 16000;                   % in Gauss/cm/s

params.os_factor = 4;               % readout oversampling amount

params.Ny_acs = 24;                 % external gre acs size ky
params.Nz_acs = 24;                 % external gre acs size kz

params.TE_acs = 6e-3;               % ACS TE
params.TR_acs = 12e-3;              % ACS TR

params.use_wave_y = use_wave_y;     % set to 1 to use sine wave on Gy 
params.use_wave_z = use_wave_z;     % set to 1 to use cosine wave on Gz
 
TEs=[3.95e-3 10.9e-3 17.85e-3 24.8e-3 31.75e-3 38.7e-3];

TR=43.5e-3;

num_echoes = length(TEs);

params.TEs = TEs;
params.TR = TR;
params.Ndummy = Ndummy;
params.N = N;
params.fov = fov;
params.flip_angle = flip_angle;
params.Tread = Tread;
params.Tpre = Tpre;
params.num_echoes = num_echoes;
params.mask = mask;
params.blip = blip;

tic
    [res_seq, res_params] = fn_generate_gre_wave_blip_shorter_label(seq, sys, params);
toc


%% Visualise sequence and output for execution

% check whether the timing of the sequence is correct
[ok, error_report]=res_seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

dirname = 'wave_gre_output/';

write_seq = 1;      % set to 1 to write .seq file


plot_start = res_params.TR * Ndummy + (Ndummy + 1) * res_params.TR;
plot_end = res_params.TR * Ndummy + (Ndummy + 5) * res_params.TR;

% visualize imaging data
res_seq.plot('TimeRange',[plot_start, plot_end])    



if write_seq
    res_seq.setDefinition('FOV', fov);
    res_seq.setDefinition('Name', 'me_gre3d');
    
    res_seq.write([dirname, 'me','_e', num2str(num_echoes),'_', num2str(1000*res(1))...
        'mmiso', '_poidisc', num2str(poidisc), '_wavey', num2str(use_wave_y), ...
        '_wavez', num2str(use_wave_z),'_blip', num2str(blip),...
        '_af',num2str(R),'_gmax',num2str(params.gZmax),'_wcycle',num2str(params.NcyclesGz), '.seq']);
end

% visualize the 3D k-space (only makes sense for low-res, otherwise one sees nothing)
if Nx<=32
    tic;
    [kfa,ta,kf,tf]=res_seq.calculateKspacePP();
    toc
    figure;plot3(kf(1,:),kf(2,:),kf(3,:));
    hold on;plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');
end


%--------------------------------------------------------------------------
%% check forbidden frequencies
%--------------------------------------------------------------------------

bay = 4;

if bay == 4
    % bay4 prisma:
    freq1 = 590;            % Hz
    bandwidth1 = 100;       % Hz

    freq2 = 1140;
    bandwidth2 = 220;
end


if bay == 2
    % bay2 terra:
    freq1 = 1100;           % Hz
    bandwidth1 = 300;       % Hz

    freq2 = 550;            % Hz
    bandwidth2 = 100;       % Hz
end


% form the Gx waveform based on the trapezoid object
GwaveX = res_params.gx.amplitude * ones(1, round(res_params.gx.flatTime / sys.gradRasterTime));

GX_rampup = linspace(0, res_params.gx.amplitude, round(res_params.gx.riseTime / sys.gradRasterTime));
GX_rampdown = linspace(res_params.gx.amplitude, 0, round(res_params.gx.fallTime / sys.gradRasterTime));

GwaveX = cat(2, GX_rampup, GwaveX, GX_rampdown);

res_params.gx.waveform = GwaveX;


h = figure(10); close(h)
figure(10), plot(res_params.gy_wave.tt, res_params.gy_wave.waveform)
figure(10), hold on, plot(res_params.gz_wave.tt, res_params.gz_wave.waveform, 'r'), axis tight
figure(10), hold on, plot(res_params.gz_wave.tt, res_params.gx.waveform, 'k'), axis tight
title('Gx, Gy, Gz waveforms during the readout')


%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------

gy_spectrum = fftc(res_params.gy_wave.waveform * use_wave_y,2);
gz_spectrum = fftc(res_params.gz_wave.waveform * use_wave_z,2);
gx_spectrum = fftc(res_params.gx.waveform,2);

points2use = 1 + length(gy_spectrum)/2 : length(gy_spectrum);   % display only one side of the spectra after the DC point

delta_freq = 1 / (length(gy_spectrum) * sys.gradRasterTime);    % in Hz



% add forbidden freq1 of the system to the plot
forbid_freq1 = zeros(1, length(points2use));

forbid_freq1_start = floor((freq1 - bandwidth1/2) / delta_freq);
forbid_freq1_end = ceil((freq1 + bandwidth1/2) / delta_freq);

% assign constant 1e7 value forbidden freq for display
forbid_freq1(1, forbid_freq1_start:forbid_freq1_end) = 1e7;     


% add forbidden freq2 of the system to the plot
forbid_freq2 = zeros(1, length(points2use));

forbid_freq2_start = floor((freq2 - bandwidth2/2) / delta_freq);
forbid_freq2_end = ceil((freq2 + bandwidth2/2) / delta_freq);

% assign constant 1e7 value forbidden freq for display
forbid_freq2(1, forbid_freq2_start:forbid_freq2_end) = 1e7;


% plot Gx, Gy, Gz and forbidden freqs together
h = figure(11); close(h)
figure(11), hold on, plot( [0:length(points2use)-1] * delta_freq, abs(gy_spectrum(points2use)), 'linewidth', 3), axis tight 
figure(11), plot( [0:length(points2use)-1] * delta_freq, abs(gz_spectrum(points2use)), 'r', 'linewidth', 3, 'linestyle', '--'), axis tight 
figure(11), plot( [0:length(points2use)-1] * delta_freq, abs(gx_spectrum(points2use)), 'k', 'linewidth', 3), axis tight 

figure(11), plot( [0:length(points2use)-1] * delta_freq, forbid_freq1, 'g'), axis tight 
figure(11), plot( [0:length(points2use)-1] * delta_freq, forbid_freq2, 'm'), axis tight,  

legend({'Gy', 'Gz', 'Gx'}, 'fontsize', 32)
ylabel('gradient spectra', 'fontsize', 32)
xlabel('Hz', 'fontsize', 32)
xlim([0,1e4])
return
