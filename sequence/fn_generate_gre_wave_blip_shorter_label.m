function [seq,params] = fn_generate_gre_wave_blip_shorter_label(seq, sys, params)

fov = params.fov;
flip_angle = params.flip_angle;

use_wave_y = params.use_wave_y;         % set to 1 to use sine wave on Gy 
use_wave_z = params.use_wave_z;         % set to 1 to use cosine wave on Gz

NcyclesGy = params.NcyclesGy;
NcyclesGz = params.NcyclesGz;

gYmax = params.gYmax;
gZmax = params.gZmax;

sYmax = params.sYmax;
sZmax = params.sZmax;

os_factor = params.os_factor;
TEs = params.TEs;
TR = params.TR;

Tread = params.Tread;
Tpre = (params.Tpre/seq.gradRasterTime)*seq.gradRasterTime;

Nx = params.N(1); 
Ny = params.N(2); 
Nz = params.N(3);              % Define FOV and resolution


TimePerSineZ = Tread / NcyclesGz;       % time per sine wave in seconds

wZ = 2*pi / TimePerSineZ;

if sZmax >= wZ * gZmax
    disp('wave amplitude is not slew limited')
    G0_Z = gZmax;
else
    disp('wave amplitude is slew limited')
    G0_Z = sZmax / wZ;
end

% check the needed rise time to reach gZmax using sZmax
gx_needed_riseTime = 2.4*G0_Z / sZmax;

disp(['needed ramp up time at least: ', num2str(1e6*gx_needed_riseTime), ' usec'])

% use this rise time for Gx gradient as a relaxed way to reach Gz max
% well within the slew limitation
gx_riseTime2use = sys.gradRasterTime * (1 + ceil(gx_needed_riseTime / sys.gradRasterTime));

disp(['ramp up time to use: ', num2str(1e6*gx_riseTime2use), ' usec'])

%--------------------------------------------------------------------------
% calculate timing: imaging data
% FOR WAVE, EXTRA GZ REWINDERS ARE IMPLEMENTED FOR MULTI-ECHO CASE TO 
% ACCOUNT FOR THE NET GZ MOMENT DUE TO COSINE WAVE
%--------------------------------------------------------------------------

num_echoes = length(TEs);
readoutSpoil = 1;

% Create non-selective pulse
[rf, rfDelay] = mr.makeBlockPulse(flip_angle*pi/180,sys,'Duration',0.2e-3);
rf_phase = 50;

% Define other gradients and ADC events
deltak=1./fov;

gx = mr.makeTrapezoid('x',sys,'FlatArea',Nx*deltak(1),'FlatTime',Tread, 'riseTime', gx_riseTime2use, 'fallTime', gx_riseTime2use); % readout gradient

adc = mr.makeAdc(Nx * os_factor,'Duration',gx.flatTime,'Delay',gx.riseTime);

adc.dwell = round(adc.dwell /seq.adcRasterTime) *seq.adcRasterTime; % ensures the timings work...



% needed only for multi-echo:
Tpre2 = ((Tpre / 3)/seq.gradRasterTime)*seq.gradRasterTime;% Tpre2 = 5*mr.calcDuration(gxPre2x);
Tpre_spoil = ((Tpre2 * 1.4)/seq.gradRasterTime)*seq.gradRasterTime;% Tpre2 = 5*mr.calcDuration(gxPre2x);

% gxPre2x = mr.makeTrapezoid('x',sys,'Area',-gx.area,'amplitude',0.8*sys.maxGrad);        % Gx prewinder -> 2x area to rewind before next echo
gxPre2x = mr.makeTrapezoid('x',sys,'Area',-gx.area,'Duration',Tpre2);        % Gx prewinder -> 2x area to rewind before next echo


gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2,'Duration',Tpre2);

try
    gxSpoil = mr.makeTrapezoid('x',sys,'Area',gx.area*(1+readoutSpoil),'Duration',Tpre2); % arbitrarilly ask the spoiler to have half the duration of the readout
catch
    gxSpoil = mr.makeTrapezoid('x',sys,'Area',gx.area*(1+readoutSpoil),'Duration',Tpre_spoil); % arbitrarilly ask the spoiler to have half the duration of the readout
end

areaY = ((0:Ny-1)-Ny/2)*deltak(2);
areaZ = ((0:Nz-1)-Nz/2)*deltak(3);

delayTEs = zeros(num_echoes,1);

% time from center of RF pulse to the start of readout? -- so this is how
% much dead time there is before readout starts
delayTEs(1) = ceil( (TEs(1) - mr.calcDuration(rf) + mr.calcRfCenter(rf) + rf.delay - mr.calcDuration(gxPre)  ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;

disp(['echo ', num2str(1), ' delay ', num2str(1e3*delayTEs(1)), ' ms'])

for t = 2:num_echoes
    delayTEs(t) = ceil( (TEs(t) - TEs(t-1) - mr.calcDuration(gx) - mr.calcDuration(gxPre2x)) / seq.gradRasterTime) * seq.gradRasterTime;

    if delayTEs(t) < 0
        disp(['echo ', num2str(t), ' cannot be fit'])
    else
        disp(['echo ', num2str(t), ' delay ', num2str(1e3*delayTEs(t)), ' ms'])
    end
end


delayTR = ceil((TR - mr.calcDuration(rf) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - mr.calcDuration(gxSpoil) - delayTEs(1))/seq.gradRasterTime)*seq.gradRasterTime;

% delay to add after the last readout and spoiler is done until end of TR

for t = 2:num_echoes
    temp = ceil((mr.calcDuration(gxPre2x) + mr.calcDuration(gx) + delayTEs(t))/seq.gradRasterTime)*seq.gradRasterTime;

    delayTR = delayTR - temp;

    if (delayTR < 0)
        disp(['delay TR: ', num2str(delayTR*1e3), ' ms'])
        msg = strcat('Please increase TR for the ', num2str(t), 'th echo');
        error(msg)
    end
end

disp(['delay TR: ', num2str(delayTR*1e3), ' ms'])

params.TEs = TEs;
params.TR = TR;

dTR=mr.makeDelay(delayTR);
%% Wave
%--------------------------------------------------------------------------
% make Gy wave gradient: sine
% consider oversampling  
%--------------------------------------------------------------------------

% Tread is equal to gx.flatTime. add delay so that Gy wave starts at the
% beginning of Gx flat top

wavepoints = round(Tread / sys.gradRasterTime);     % number of readout points
T_wavepoints = sys.gradRasterTime;                  % grad raster time in seconds

TimePerSineY = (Tread / sys.gradRasterTime) * T_wavepoints / NcyclesGy;     % time per sine wave in seconds

wY = 2*pi / TimePerSineY;

if sYmax >= wY * gYmax
    G0_Y = gYmax;

    disp(['wave amplitude is not slew limited, using gY = ', num2str(G0_Y*10), ' mT/m'])
else
    
    G0_Y = sYmax / wY;
    
    disp(['wave amplitude is slew limited, using gY = ', num2str(G0_Y*10), ' mT/m'])
end

SO_Y = G0_Y * wY;

% add one last point to go back to zero Gy (wavepoints+1 total)
% this seems to be done in the siemens code as well
GradTimePointsY = [0:wavepoints] * T_wavepoints;        

scaling_factor = sys.gamma * 1e-2;                          % convert from G/cm to T/m by multiplying by 1e-2, then to Hz/m by multiplying by gamma

GwaveY = G0_Y * sin(wY * GradTimePointsY) * scaling_factor; % Gy amplitude in Hz/m


disp(['Gx amplitude: ', num2str(max(abs(gx.amplitude)*1e-3)), ' kHz/m'])
disp(['Gy amplitude: ', num2str(max(abs(GwaveY)*1e-3)), ' kHz/m'])


% pad Gy waveform at the beginning and end with zeroes to account for rise and fall time of trapezoid
num_pad_pre = round(gx_riseTime2use / sys.gradRasterTime);


% -1 to account for added point at the end to make sure sine wave goes back
% to zero due to [0:wavepoints] instead of [0:wavepoints-1]
num_pad_post = round(gx_riseTime2use / sys.gradRasterTime) - 1;     

GwaveY = cat(2, zeros(1,num_pad_pre), GwaveY, zeros(1,num_pad_post));

disp(['num points Gy:', num2str(length(GwaveY))])


% form the Gx waveform for display
GwaveX = gx.amplitude * ones(1, round(gx.flatTime / sys.gradRasterTime));

GX_rampup = linspace(0, gx.amplitude, round(gx.riseTime / sys.gradRasterTime));
GX_rampdown = linspace(gx.amplitude, 0, round(gx.fallTime / sys.gradRasterTime));

GwaveX = cat(2, GX_rampup, GwaveX, GX_rampdown);

disp(['num points Gx:', num2str(length(GwaveX))])


% plot Gx and Gy gradients together
h = figure(20), close(h);
figure(20), plot(GwaveY,'o'), xlabel('sec'), ylabel('kHz/m')
figure(20), hold on, plot(GwaveX, 'rx'), xlabel('sec'), ylabel('kHz/m')


% form the Gy object
gy_wave = mr.makeArbitraryGrad('y', GwaveY, sys);
gy_wave.id = seq.registerGradEvent(gy_wave);
gy_wave.first = 0;
gy_wave.last = 0;

% gy_balance_wav = mr.makeTrapezoid('y',sys,'Area',-gy_wave.area,'Duration',Tpre);


%--------------------------------------------------------------------------
% make Gz wave gradient: cosine
% consider oversampling  
%--------------------------------------------------------------------------


TimePerSineZ = (Tread / sys.gradRasterTime) * T_wavepoints / NcyclesGz;     % time per sine wave in seconds

wZ = 2*pi / TimePerSineZ;

if sZmax >= wZ * gZmax
    G0_Z = gZmax;

    disp(['wave amplitude is not slew limited, using gZ = ', num2str(G0_Z * 10), ' mT/m'])
else
    
    G0_Z = sZmax / wZ;

    disp(['wave amplitude is slew limited, using gZ = ', num2str(G0_Z * 10), ' mT/m'])
end

SO_Z = G0_Z * wZ;


GradTimePointsZ = [0:wavepoints-1] * T_wavepoints;       % add one last point to go back to zero Gy (wavepoints+1 total)

scaling_factor = sys.gamma * 1e-2;  % convert from G/cm to T/m by multiplying by 1e-2, then to Hz/m by multiplying by gamma


GwaveZ_flat = G0_Z * cos(wZ * GradTimePointsZ) * scaling_factor;

GZ_rampup = linspace(0, G0_Z * scaling_factor, round(gx.riseTime / sys.gradRasterTime));
GZ_rampdown = linspace(G0_Z * scaling_factor, 0, round(gx.fallTime / sys.gradRasterTime));

GwaveZ = cat(2, GZ_rampup, GwaveZ_flat, GZ_rampdown);

disp(['num points Gx:', num2str(length(GwaveZ))])


% plot Gx and Gy gradients together
h = figure(20), close(h);
figure(20), plot( GwaveY, 'linewidth', 3), xlabel('sec'), ylabel('kHz/m')
figure(20), hold on, plot( GwaveX, 'r'), xlabel('sec'), ylabel('kHz/m')
figure(20), hold on, plot( GwaveZ, 'k', 'linewidth', 3), xlabel('sec'), ylabel('kHz/m'), 


% form the Gz object
gz_wave = mr.makeArbitraryGrad('z', GwaveZ,sys);

% manually set the first and the end to be zero
gz_wave.first = 0;
gz_wave.last = 0;
gz_wave.id = seq.registerGradEvent(gz_wave);


gz_balance_wav = mr.makeTrapezoid('z',sys,'Area', -gz_wave.area,'Duration',Tpre2);

%% Other parameters

Tpre = params.Tpre;
Ndummy = params.Ndummy;

if params.blip
    blips = params.blip_shift;
end

% Make trapezoids for inner loop to save computation
clear gyPre gyReph;

for iY=1:Ny
    gyPre(iY) = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',Tpre2);
    gyReph(iY) = mr.makeTrapezoid('y','Area',-areaY(iY),'Duration',Tpre2);
end

for iZ=1:Nz
    gzPre(iZ) = mr.makeTrapezoid('z','Area',areaZ(iZ),'Duration',Tpre2);
    gzReph(iZ) = mr.makeTrapezoid('z','Area',-areaZ(iZ),'Duration',Tpre2);
end


% combine blip and gz balancing gradnets
if params.blip
    for u = 2:num_echoes
        gz_combined(u) = mr.makeTrapezoid('z',sys,'Area', -gz_wave.area+gzPre(Nz/2+1+round(blips(2,u-1))).area,'Duration',Tpre2);
    end
end

% Drive magnetization to the steady state
for iY=1:Ndummy
    % RF
    rf.phaseOffset = mod(rf_phase*(iY^2+iY+2)*pi/180,2*pi);
    seq.addBlock(rf,rfDelay);
    % Gradients
    seq.addBlock(gxPre,gyPre(round(Ny/2)),gzPre(round(Nz/2)));
    seq.addBlock(mr.makeDelay(delayTEs(1)));

    if use_wave_y && use_wave_z
        seq.addBlock(gx, gy_wave, gz_wave);                   
    end
        
    if use_wave_y && ~use_wave_z
        seq.addBlock(gx, gy_wave);                   
    end
    
    if ~use_wave_y && use_wave_z
        seq.addBlock(gx, gz_wave);
    end
    
    
    if ~use_wave_y && ~use_wave_z
        seq.addBlock(gx);
    end

    for t = 2:num_echoes
        seq.addBlock(mr.makeDelay(delayTEs(t)));    % delay until readout
        
        if ~use_wave_y && use_wave_z
           seq.addBlock(gxPre2x, gz_balance_wav);
        end
                
        if ~use_wave_z && use_wave_y
           seq.addBlock(gxPre2x);
        end
                
        if use_wave_y && use_wave_z
            seq.addBlock(gxPre2x, gz_balance_wav);                   
        end
                   
        if ~use_wave_y && ~use_wave_z
            seq.addBlock(gxPre2x);
        end
        
        if use_wave_y && use_wave_z
            seq.addBlock(gx, gy_wave, gz_wave);                   
        end
        
        if use_wave_y && ~use_wave_z
            seq.addBlock(gx, gy_wave);                   
        end
    
        if ~use_wave_y && use_wave_z
            seq.addBlock(gx, gz_wave); 
        end
    
        if ~use_wave_y && ~use_wave_z
            seq.addBlock(gx);
        end
    end
    
    if use_wave_z
        areaz2 = -areaZ(round(Nz/2))-gz_wave.area;
    else
        areaz2 = -areaZ(round(Nz/2));
    end
    
    areay2 = -areaY(round(Ny/2));
    
    gyReph2 = mr.makeTrapezoid('y','Area',areay2,'Duration',Tpre2);
    gzReph2 = mr.makeTrapezoid('z','Area',areaz2,'Duration',Tpre2);
    
    seq.addBlock(gyReph2, gzReph2, gxSpoil);
    seq.addBlock(dTR);
end

mask = params.mask;
shift_scale_y = 1;
shift_scale_z = 1;


% Loop over phase encodes and define sequence blocks
for iZ=1:Nz
    for iY=1:Ny

        if(abs(mask(iY,iZ)) > 5e-1)
            % RF spoiling
            rf.phaseOffset = mod(rf_phase*(iY^2+iY+2)*pi/180,2*pi);
            adc.phaseOffset = rf.phaseOffset;
        
            % Excitation
            seq.addBlock(rf,rfDelay);
            
            blip_y = zeros(num_echoes-1,1);
            blip_z = zeros(num_echoes-1,1);
            
            sum_blip_z = 0;
            sum_blip_y = 0;
        
            % Encoding
            % add labels
            PrePhaserBlockContents={gxPre,gyPre(iY),gzPre(iZ)};
            
            PrePhaserBlockContents= { PrePhaserBlockContents{:} , mr.makeLabel('SET','SET', 0) }; % set the echo counters
            
            PrePhaserBlockContents= { PrePhaserBlockContents{:}, mr.makeLabel('SET','LIN',iY-1), mr.makeLabel('SET','PAR', iZ-1) }; % set the lin and partitions
            
            seq.addBlock(PrePhaserBlockContents{:});

            seq.addBlock(mr.makeDelay(delayTEs(1)));

            if use_wave_y && use_wave_z
                seq.addBlock(gx, gy_wave, gz_wave, adc);                   
            end

            if use_wave_y && ~use_wave_z
                seq.addBlock(gx, gy_wave, adc);                   
            end

            if ~use_wave_y && use_wave_z
                seq.addBlock(gx, gz_wave, adc);                 
            end

            if ~use_wave_y && ~use_wave_z
                seq.addBlock(gx, adc);                 
            end

            for t = 2:num_echoes
                seq.addBlock(mr.makeDelay(delayTEs(t)),mr.makeLabel('SET','SET', t-1));    % delay until readout
                                
                % compensation gradients for Wave readout (only z direction needs to be compensated)
                if use_wave_z
                    if params.blip
                        if (abs(blips(1,t-1)) < 1e-13) && (abs(blips(2,t-1)) > 1e-2)
                            
                            blip_z(t-1) = round(blips(2,t-1)); 
                            sum_blip_z = sum_blip_z + blip_z(t-1);
                            
                            if (((iZ + sum_blip_z) > Nz) || ((iZ + sum_blip_z) < 1))
                                blip_z(t-1) = 0;
                                RephaserBlockContents = {gxPre2x,gz_balance_wav};
                                sum_blip_z = sum(blip_z(1:(t-1)),1);
                            else
                                RephaserBlockContents = {gxPre2x,gz_combined(t)};
                            end
                                                        
                            RephaserBlockContents = {RephaserBlockContents{:}, ...
                            mr.makeLabel('SET','LIN', iY-1 + round(sum_blip_y)), ...
                            mr.makeLabel('SET','PAR', iZ-1 + round(sum_blip_z))};
                            seq.addBlock(RephaserBlockContents{:});
                        end
                
                        if (abs(blips(1,t-1)) > 1e-2) && (abs(blips(2,t-1)) < 1e-13)
                            
                            blip_y(t-1) = round(blips(1,t-1));
                            sum_blip_y = sum_blip_y + blip_y(t-1);
                            
                            if (((iY +sum_blip_y) > Ny) || ((iY + sum_blip_y) < 1))
                                blip_y(t-1) = 0;
                                RephaserBlockContents = {gxPre2x,gz_balance_wav};
                                sum_blip_y = sum(blip_y(1:(t-1)),1);
                            else
                                RephaserBlockContents = {gxPre2x,gyPre(Ny/2+1+round(blip_y(t-1))), gz_balance_wav};
                            end
                            
                            RephaserBlockContents = {RephaserBlockContents{:}, ...
                            mr.makeLabel('SET','LIN', iY-1 + round(sum_blip_y)), ...
                            mr.makeLabel('SET','PAR', iZ-1 + round(sum_blip_z))};
                            seq.addBlock(RephaserBlockContents{:});
                        end
                
                        if (abs(blips(1,t-1)) > 1e-2) && abs((blips(2,t-1)) > 1e-2)

                            blip_y(t-1) = round(blips(1,t-1));
                            blip_z(t-1) = round(blips(2,t-1));
                            
                            sum_blip_z = sum_blip_z + blip_z(t-1);
                            sum_blip_y = sum_blip_y + blip_y(t-1);
                            
                            if (((iZ + sum_blip_z) > Nz) || ((iZ + sum_blip_z) < 1))
                                blip_z(t-1) = 0;
                                sum_blip_z = sum(blip_z(1:(t-1)),1);
                            end
                            
                            if (((iY + sum_blip_y) > Ny) || ((iY + sum_blip_y) < 1))
                                blip_y(t-1) = 0;
                                sum_blip_y = sum(blip_y(1:(t-1)),1);
                            end
                            
                            if ~blip_z(t-1) && ~blip_y(t-1)
                                RephaserBlockContents = {gxPre2x,gz_balance_wav};
                            end
                            
                            if blip_z(t-1) && ~blip_y(t-1)
                                RephaserBlockContents = {gxPre2x,gz_combined(t)};
                            end
                            
                            if ~blip_z(t-1) && blip_y(t-1)
                                RephaserBlockContents = {gxPre2x,gyPre(Ny/2+1+round(blip_y(t-1))), gz_balance_wav};
                            end
                            
                            if blip_z(t-1) && blip_y(t-1)
                                RephaserBlockContents = {gxPre2x,gyPre(Ny/2+1+round(blip_y(t-1))), gz_combined(t)};
                            end

                            RephaserBlockContents = {RephaserBlockContents{:}, ...
                            mr.makeLabel('SET','LIN', iY-1 + round(sum_blip_y)), ...
                            mr.makeLabel('SET','PAR', iZ-1 + round(sum_blip_z))};
                            seq.addBlock(RephaserBlockContents{:});
                        end
                        
                    else
                        
                        RephaserBlockContents = {gxPre2x,gz_balance_wav};
                        RephaserBlockContents = {RephaserBlockContents{:}, ...
                        mr.makeLabel('SET','LIN', iY-1), ...
                        mr.makeLabel('SET','PAR', iZ-1)};
                        seq.addBlock(RephaserBlockContents{:});
                    end 
                end
                
                if ~use_wave_z
                   if params.blip                       
                        if (abs(blips(1,t-1)) < 1e-13) && (abs(blips(2,t-1)) > 1e-2)
                            blip_z(t-1) = round(blips(2,t-1));
                            
                            sum_blip_z = sum_blip_z + blip_z(t-1);
                            
                            if (((iZ + sum_blip_z) > Nz) || ((iZ + sum_blip_z) < 1))
                                blip_z(t-1) = 0;
                                RephaserBlockContents = {gxPre2x};
                                sum_blip_z = sum(blip_z(1:(t-1)),1);
                            else
                                RephaserBlockContents = {gxPre2x,gzPre(Nz/2+1+round(blip_z(t-1)))};
                            end
                                      
                            RephaserBlockContents = {RephaserBlockContents{:}, ...
                            mr.makeLabel('SET','LIN', iY-1 + round(sum_blip_y)), ...
                            mr.makeLabel('SET','PAR', iZ-1 + round(sum_blip_z))};
                            seq.addBlock(RephaserBlockContents{:});
                        end
                
                        if (abs(blips(1,t-1)) > 1e-2) && (abs(blips(2,t-1)) < 1e-13)
                            blip_y(t-1) = round(blips(1,t-1));
                            
                            sum_blip_y = sum_blip_y + blip_y(t-1);
                            
                            if (((iY + sum_blip_y) > Ny) || ((iY + sum_blip_y) < 1))
                                blip_y(t-1) = 0;
                                RephaserBlockContents = {gxPre2x};
                                sum_blip_y = sum(blip_y(1:(t-1)),1);
                            else
                                RephaserBlockContents = {gxPre2x,gyPre(Ny/2+1+round(blip_y(t-1)))};
                            end

                            RephaserBlockContents = {RephaserBlockContents{:}, ...
                            mr.makeLabel('SET','LIN', iY-1 + round(sum_blip_y)), ...
                            mr.makeLabel('SET','PAR', iZ-1 + round(sum_blip_z))};
                            seq.addBlock(RephaserBlockContents{:});
                        end
                
                        if (abs(blips(1,t-1)) > 1e-2) && abs((blips(2,t-1)) > 1e-2)

                            blip_y(t-1) = round(blips(1,t-1));
                            blip_z(t-1) = round(blips(2,t-1));
                            
                            sum_blip_y = sum_blip_y + blip_y(t-1);
                            sum_blip_z = sum_blip_z + blip_z(t-1);
                            
                            if (((iY + sum_blip_y) > Ny) || ((iY + sum_blip_y) < 1))
                                blip_y(t-1) = 0;
                                sum_blip_y = sum(blip_y(1:(t-1)),1);
                            end
                            
                            if (((iZ + sum_blip_z) > Nz) || ((iZ + sum_blip_z) < 1))
                                blip_z(t-1) = 0;
                                sum_blip_z = sum(blip_z(1:(t-1)),1);
                            end

                            if ~blip_z(t-1) && ~blip_y(t-1)
                                RephaserBlockContents = {gxPre2x};
                            end
                            
                            if blip_z(t-1) && ~blip_y(t-1)
                                RephaserBlockContents = {gxPre2x,gzPre(Nz/2+1+round(blip_z(t-1)))};
                            end
                            
                            if ~blip_z(t-1) && blip_y(t-1)
                                RephaserBlockContents = {gxPre2x,gyPre(Ny/2+1+round(blip_y(t-1)))};
                            end
                            
                            if blip_z(t-1) && blip_y(t-1)
                                RephaserBlockContents = {gxPre2x,gyPre(Ny/2+1+round(blip_y(t-1))),gzPre(Nz/2+1+round(blip_z(t-1)))};
                            end
                            
                            RephaserBlockContents = {RephaserBlockContents{:}, ...
                            mr.makeLabel('SET','LIN', iY-1 + round(sum_blip_y)), ...
                            mr.makeLabel('SET','PAR', iZ-1 + round(sum_blip_z))};
                            seq.addBlock(RephaserBlockContents{:});
                        end
                    else
                        RephaserBlockContents = {gxPre2x};
                        RephaserBlockContents = {RephaserBlockContents{:}, ...
                        mr.makeLabel('SET','LIN', iY-1), ...
                        mr.makeLabel('SET','PAR', iZ-1)};
                        seq.addBlock(RephaserBlockContents{:});
                   end

                end

                
                % Wave readouts for later echos
                if use_wave_y && use_wave_z
                    seq.addBlock(gx, gy_wave, gz_wave, adc);                   
                end

                if use_wave_y && ~use_wave_z
                    seq.addBlock(gx, gy_wave, adc);
                end

                if ~use_wave_y && use_wave_z
                    seq.addBlock(gx, gz_wave, adc); 
                end

                if ~use_wave_y && ~use_wave_z
                    seq.addBlock(gx,adc);
                end
            end
            
            if use_wave_z
                areaz2 = -areaZ(iZ)-gz_wave.area;
            else
                areaz2 = -areaZ(iZ);
            end
            
            areay2 = -areaY(iY);
    
            if params.blip
                areay2 = areay2 - areaY(Ny/2+1+round(sum(blip_y(1:num_echoes-1),1)));
                areaz2 = areaz2 - areaZ(Nz/2+1+round(sum(blip_z(1:num_echoes-1),1)));
            end

            
            gyReph2 = mr.makeTrapezoid('y','Area',areay2,'Duration',Tpre2);
            gzReph2 = mr.makeTrapezoid('z','Area',areaz2,'Duration',Tpre2);
            % 
            seq.addBlock(gyReph2, gzReph2, gxSpoil);
            seq.addBlock(dTR);

        end
    end
end

params.mask = mask;
params.gx = gx;
params.gy_wave = gy_wave;
params.gz_wave = gz_wave;

fprintf('Sequence ready\n');
