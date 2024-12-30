% This present the MATLAB codes for LPI radar signal generation and Choi-Williams distribution-based time-frequency image generation.

% The codes of Choi-Williams distribution have been achieved from the
% Time-Frequency Toolbox in http://tftb.nongnu.org/
% The TFTB is distributed under the terms of the GNU Public Licence.

% % Note: before running this codes, please extract the meterials.zip file into the same folder with this code.
% addpath(genpath("\tftb-0.2")); % Modify the pathname in your pc
% addpath 'waveform-types'

%% initial parameters configurations
fs = 100e6; % sample frequency
A = 1;      % amplitude
waveforms = {'Rect','LFM','Costas','Barker','Frank','P1','P2','P3','P4','T1','T2','T3','T4','ClutterOnly'};   % 13 LPI waveform codes
datasetCWD = 'dataset-CWD-64';

for i = 1 : length(waveforms)
    % create the folders for dataset storage
    mkdir(fullfile(datasetCWD,waveforms{i}));
end

fps = 1000;  % the number of signal per SNR per waveform codes
% filter configuration
g=kaiser(63,0.5);
h=kaiser(63,0.5);
imgSize = 64;

%% initial channels
atenuationset = -20:2:0; % path gain range
delayset = (1:50:1000)*1e-9; % path delay range
atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
delays = [0,delayset(randperm(length(delayset),5))]; % path delays
multipathChannel = comm.RayleighChannel(...
    'SampleRate', fs, ...
    'PathDelays', delays, ...
    'AveragePathGains', atenuations, ...
    'MaximumDopplerShift', randi([5 400]));

SNR = 1 : 1 : 10;     % snr range

for n = 1 : length(SNR)
    disp(['SNR = ',sprintf('%+02d',SNR(n))])
snr_num=SNR(n);
    for k = 1 : length(waveforms)
        waveform = waveforms{k};
        switch waveform
            case 'Rect'
                disp(['Generating ',waveform, ' waveform ...']);
                fc = linspace(fs/6,fs/5,fps);
                fc = fc(randperm(fps));
                N = linspace(512,1920,fps);
                N = round(N(randperm(fps)));
                waveformfolderCWD = fullfile(datasetCWD,waveform);

                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);   % Make a copy

                    release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % Uniformly random path gains generation
                    delays = [0,delayset(randperm(length(delayset),5))]; % Uniformly random path delays generation
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    wav = type_Rect(N(idx),fs,A,fc(idx));
                    wav = localChannel(wav');       % passing over multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('rect_snr_%02d_no%05d.png',snr_num,idx)));

                end
            case 'LFM'
                disp(['Generating ',waveform, ' waveform ...']);
                fc = linspace(fs/6,fs/5,fps);
                fc=fc(randperm(fps));
                B = linspace(fs/20, fs/16, fps);
                B = B(randperm(fps));
                N = linspace(512,1920,fps);
                N=round(N(randperm(fps)));
                sweepDirections = {'Up','Down'};
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    wav = type_LFM(N(idx),fs,A,fc(idx),B(idx),sweepDirections{randi(2)});
                    wav = localChannel(wav'); % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('lfm_snr_%02d_no%05d.png',snr_num,idx)));

                end
            case 'Costas'
                disp(['Generating ',waveform, ' waveform ...']);
                Lc = [3, 4, 5, 6];
                fcmin = linspace(fs/30,fs/24,fps);
                fcmin=fcmin(randperm(fps));
                N = linspace(512,1920,fps);
                N=round(N(randperm(fps)));
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    NumHop = randperm(Lc(randi(4)));
                    wav = type_Costas(N(idx), fs, A, fcmin(idx), NumHop);
                    wav = localChannel(wav'); % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('costas_snr_%02d_no%05d.png',snr_num,idx)));

                end

            case 'Barker'
                disp(['Generating ',waveform, ' waveform ...']);
                Lc = [7,11,13];
                fc = linspace(fs/6,fs/5,fps);
                fc=fc(randperm(fps));
                Ncc = 20:24;
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    Bar = Lc(randi(3));
                    if Bar == 7
                        phaseCode = [0 0 0 1 1 0 1]*pi;
                    elseif Bar == 11
                        phaseCode = [0 0 0 1 1 1 0 1 1 0 1]*pi;
                    elseif Bar == 13
                        phaseCode = [0 0 0 0 0 1 1 0 0 1 0 1 0]*pi;
                    end
                    wav = type_Barker(Ncc, fs, A, fc(idx), phaseCode);
                    wav = localChannel(wav'); % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('barker_snr_%02d_no%05d.png',snr_num,idx)));

                end
            case 'Frank'
                disp(['Generating ',waveform, ' waveform ...']);
                fc = linspace(fs/6,fs/5,fps);
                fc=fc(randperm(fps));
                Ncc = [3,4,5];
                M = [6, 7, 8];
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    wav = type_Frank(Ncc(randi(3)), fs, A, fc(idx), M(randi(3)));
                    wav = localChannel(wav'); % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('frank_snr_%02d_no%05d.png',snr_num,idx)));

                end
            case 'P1'
                disp(['Generating ',waveform, ' waveform ...']);
                fc = linspace(fs/6,fs/5,fps);
                fc=fc(randperm(fps));
                Ncc = [3,4,5];
                M = [6, 7, 8];
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    wav = type_P1(Ncc(randi(3)), fs, A, fc(idx), M(randi(3)));
                    wav = localChannel(wav'); % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('p1_snr_%02d_no%05d.png',snr_num,idx)));

                end
            case 'P2'
                disp(['Generating ',waveform, ' waveform ...']);
                fc = linspace(fs/6,fs/5,fps);
                fc=fc(randperm(fps));
                Ncc = [3,4,5];
                M = [6, 8];
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    wav = type_P2(Ncc(randi(3)), fs, A, fc(idx), M(randi(2)));
                    wav = localChannel(wav'); % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('p2_snr_%02d_no%05d.png',snr_num,idx)));

                end
            case 'P3'
                disp(['Generating ',waveform, ' waveform ...']);
                fc = linspace(fs/6,fs/5,fps);
                fc=fc(randperm(fps));
                Ncc = [3,4,5];
                p = [36, 49, 64];
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    wav = type_P3(Ncc(randi(3)), fs, A, fc(idx), p(randi(3)));
                    wav = localChannel(wav'); % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('p3_snr_%02d_no%05d.png',snr_num,idx)));

                end

            case 'P4'
                disp(['Generating ',waveform, ' waveform ...']);
                fc = linspace(fs/6,fs/5,fps);
                fc=fc(randperm(fps));
                Ncc = [3,4,5];
                p = [36, 49, 64];
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    wav = type_P4(Ncc(randi(3)), fs, A, fc(idx), p(randi(3)));
                    wav = localChannel(wav'); % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('p4_snr_%02d_no%05d.png',snr_num,idx)));

                end
            case 'T1'
                disp(['Generating ',waveform, ' waveform ...']);
                fc = linspace(fs/6,fs/5,fps);
                fc=fc(randperm(fps));
                Ng = [4,5,6];
                N = linspace(512,1920,fps);
                N=round(N(randperm(fps)));
                Nps = 2;
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    wav = type_T1(fs, A, fc(idx),Nps,Ng(randi(3)));
                    wav = localChannel(wav'); % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('t1_snr_%02d_no%05d.png',snr_num,idx)));

                end

            case 'T2'
                disp(['Generating ',waveform, ' waveform ...']);
                fc = linspace(fs/6,fs/5,fps);
                fc=fc(randperm(fps));
                Ng = [4,5,6];
                Nps = 2;
                N = linspace(512,1920,fps);
                N=round(N(randperm(fps)));
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    wav = type_T2(fs, A, fc(idx),Nps,Ng(randi(3)));
                    wav = localChannel(wav'); % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('t2_snr_%02d_no%05d.png',snr_num,idx)));

                end

            case 'T3'
                disp(['Generating ',waveform, ' waveform ...']);
                fc = linspace(fs/6,fs/5,fps);
                fc=fc(randperm(fps));
                B = linspace(fs/20,fs/10,fps);
                B = B(randperm(fps));
                Ng = [4,5,6];
                N = linspace(512,1920,fps);
                N=round(N(randperm(fps)));
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    wav = type_T3(N(idx), fs, A, fc(idx), Nps,B(idx));
                    wav = localChannel(wav'); % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('t3_snr_%02d_no%05d.png',snr_num,idx)));

                end

            case 'T4'
                disp(['Generating ',waveform, ' waveform ...']);
                fc = linspace(fs/6,fs/5,fps);
                fc=fc(randperm(fps));
                B = linspace(fs/20,fs/10,fps);
                B = B(randperm(fps));
                Ng = [4,5,6];
                N = linspace(512,1920,fps);
                N=round(N(randperm(fps)));
                waveformfolderCWD = fullfile(datasetCWD,waveform);
                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0,atenuationset(randperm(length(atenuationset),5))]; % path gains
                    delays = [0,delayset(randperm(length(delayset),5))]; % path delays
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);
                    wav = type_T4(N(idx), fs, A, fc(idx), Nps,B(idx));
                    wav = localChannel(wav');       % adding multipath channels
                    wav = awgn(wav,SNR(n),'measured');
                    t=1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD,fullfile(waveformfolderCWD,sprintf('t4_snr_%02d_no%05d.png',snr_num,idx)));

                end

            case 'ClutterOnly'
                disp('Generating Clutter-Only waveforms...');
                waveformfolderCWD = fullfile(datasetCWD, 'ClutterOnly');
                fc = linspace(fs/6, fs/5, fps);  % not actually used, but to match code structure
                N = linspace(512,1920,fps);
                N = round(N(randperm(fps)));

                parfor idx = 1:fps
                    localChannel = clone(multipathChannel);release(localChannel);
                    atenuations = [0, atenuationset(randperm(length(atenuationset),5))];
                    delays = [0, delayset(randperm(length(delayset),5))];
                    localChannel.PathDelays = delays;
                    localChannel.AveragePathGains = atenuations;
                    localChannel.MaximumDopplerShift = randi([10 1000]);

                    % 1) random noise
                    wav = randn(N(idx), 1);
                    % 2) pass through channel
                    wav = localChannel(wav);
                    % 3) optionally add AWGN at SNR(n)
                    wav = awgn(wav, SNR(n), 'measured');

                    t = 1:length(wav);
                    [CWD_TFD,~,~] = FTCWD(wav,t,1024,g,h,1,0,imgSize);
                    imwrite(CWD_TFD, ...
                        fullfile(waveformfolderCWD, ...
                        sprintf('clutter-snr%02d-no%05d.png', n, idx)));
                end

            otherwise
                disp('Done!')
        end

    end

end


