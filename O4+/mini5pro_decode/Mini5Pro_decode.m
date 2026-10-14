%% Mini5Pro V4.2 Syntax-Clean + FAST + Live Constellations (pure script)
% Resolve V4.1 MATLAB parser issues by normalizing continuation lines and whitespace.
% Retain Mini5Pro V4.0 ZC, pilots, CP/CFO, QAM/EVM processing, and file output logic.
% Address slow full-record high-order FIR filtering and 17 repeated sliding-energy calculations;
% also improve late plotting and figure lifetime without changing ZC root, OFDM bins, or 13-symbol structure.
% Date: 2026-10-14
% 2026-09-22
%
% Purpose:
%   raw MAT IQ
% Locate the Mini5Pro/O4 1201-tone video-transmission channel center automatically.
% DDC and decimation to 30.72 MS/s.
% Fractional CFO and integer-bin CFO correction.
% ZC (L=1201, root=232) burst detection.
% Structured timing of 27 symbols, CP-based CFO, and fine CFO from ZC.
% Time interpolation of ZC-based channel estimates.
% m1 QPSK and geometric 16/64QAM decisions for other data symbols.
% Classify 13-symbol blocks as NORMAL_64, ALL_QPSK, or MIXED.
% EVM, CP-correlation, and m1 recurrence validation.
% Repeated 64QAM hard-level consistency for equal rel positions in the same file.
% Single/multi-file constellation figures, CSV and MAT outputs.
%
% This pure-script directory batch processor follows the previously evaluated mini5pro_rx_v3.py / mini5pro_rxlib.py receiver.
% Do not restore the old v2 per-resource-element phase DP or claim geometric hard levels are verified DJI bit labels.
%
% High-confidence parameters from the previous Mini5Pro v3 receiver and a 46-file evaluation:
%   Fs_out           = 30.72 MHz
%   NFFT             = 2048
%   Delta-f          = 15 kHz
% active K         = -600:600 (1201 RE including DC)
%   ZC               = L=1201, root=232
% ZC/burst pattern = five ZC symbols; approximate peak spacings [6848 15360 13168 15360]
%   burst rel        = 0..26
%   ZC rel           = [0 3 10 16 23]
% data rel         = remaining 22 positions
% long CP rel      = [1 8 15 22] use CP=160; other data positions use CP=144
%   base step        = NFFT+144 = 2192 samples
% m1 QPSK          = rel=1 / rel=14; P131/P95 recurrence checks, respectively
% 46-file result   = 7 ALL_QPSK blocks among 272; all seven occurred in Block 1
%
% Not yet established:
% the actual DJI 64QAM six-bit labeling
%   - descrambler / interleaver / FEC
% Accordingly, the primary outputs are complex soft symbols and geometric I/Q levels.
%
% Dependency: MATLAB Signal Processing Toolbox (fir1, kaiser, etc.).
% This is a script and does not define MATLAB functions.
clc;
clear;
close all;
%% ================================================================
% 0. User configuration
% ================================================================
% Provide a MAT file or directory; directory mode searches for *.mat recursively.
inputPath = 'E:\DJI_figure\data\data\O4+\mini5pro\mini5pro_飞行_O4+\iq200ms_20260128_235329_430_000186.mat';
% Output directory; when empty, create Mini5Pro_V4_output beside the input.
outputRoot = '';
% Draw constellations for each file and the combined set of files.
PLOT_CONSTELLATION = true;
% Immediately display a lightweight raw-IQ time-frequency overview (not full STFT).
PLOT_OVERVIEW = true;
% Maximum constellation points displayed; processing and statistics still use all points.
MAX_CONSTELLATION_POINTS = 30000;
% Live constellation refresh interval in bursts (always refresh the first and last).
CONSTELLATION_UPDATE_EVERY = 3;
% Save each file's recovered Z/H/levels/modulation data to MAT.
SAVE_PER_FILE_MAT = true;
% Save 1200 non-DC I/Q hard levels per 64QAM symbol to CSV.
% May produce large CSV files; disabled by default, enable for bitstream analysis.
SAVE_HARD_LEVEL_CSV = false;
% Front-end and ZC detection thresholds.
ZC_RHO_MIN          = 0.60;
INTEGER_ZC_RHO_MIN  = 0.80;
EDGE_CONTRAST_MINDB = 10.0;
% Warn in the console if m1 recurrence falls below this threshold.
M1_WARN = 0.995;
%% ================================================================
% 1. Input file list and output directory
% ================================================================
if isempty(outputRoot)
    if isfolder(inputPath)
        outputRoot = fullfile(inputPath,'Mini5Pro_V4_output');
    else
        outputRoot = fullfile(fileparts(inputPath),'Mini5Pro_V4_output');
    end
end
if ~exist(outputRoot,'dir'), mkdir(outputRoot); end
if isfolder(inputPath)
    dd = dir(fullfile(inputPath,'**','*.mat'));
    matFiles = {};
    for ii=1:numel(dd)
        ff = fullfile(dd(ii).folder,dd(ii).name);
        % Prevent re-reading MAT files created by this script in future directory runs.
        if startsWith(ff,outputRoot,'IgnoreCase',true)
            continue;
        end
        matFiles{end+1,1} = ff; %#ok<SAGROW>
    end
else
    matFiles = {inputPath};
end
if isempty(matFiles)
    error('No MAT files found: %s',inputPath);
end
%% ================================================================
% 2. Mini5Pro PHY constants
% ================================================================
FS_OUT = 30.72e6;
NFFT   = 2048;
SCS    = FS_OUT/NFFT;             % 15 kHz
L      = 1201;
ROOT   = 232;
HALF   = 600;
K      = (-HALF:HALF).';
KIDX   = mod(K,NFFT)+1;           % unshifted FFT bin indices
NODC   = (K~=0);
ANCHOR_REL = [0 3 10 16 23];     % All ZC anchor positions.
CFO_ZC_REL = [3 10 16 23];       % Four ZC positions used for fine CFO estimation and coherence screening.
LONG_CP_REL = [1 8 15 22];
CP_N = 144;
CP_L = 160;
% m1 recurrence taps
TAPS131 = [29 30 56 57 58 60 61 62 69 98 99 125 126 127 129 130 131];
TAPS95  = [29 30 33 56 57 58 60 61 63 89 90 91 93 94 95];
% ZC frequency/time template
nzc = (0:L-1).';
ZC = exp(-1j*pi*ROOT*nzc.*(nzc+1)/L);
Gzc = complex(zeros(NFFT,1));
Gzc(KIDX) = ZC;
tZ = ifft(Gzc);
tZ = tZ/sqrt(sum(abs(tZ).^2)+eps);
% ideal constellations: geometry only
lev64 = (-7:2:7).'/sqrt(42);
[II64,QQ64] = ndgrid(lev64,lev64);
C64 = II64(:)+1j*QQ64(:);
lev16 = [-3 -1 1 3].'/sqrt(10);
[II16,QQ16] = ndgrid(lev16,lev16);
C16 = II16(:)+1j*QQ16(:);
CQ = [1+1j;1-1j;-1+1j;-1-1j]/sqrt(2);
%% ================================================================
% 3. Global result containers
% ================================================================
% file-level
F_file={}; F_state={}; F_fs=[]; F_fc=[]; F_nsamp=[]; F_duration=[];
F_center_rel_hz=[]; F_center_abs_hz=[]; F_edge_db=[]; F_nburst=[];
% burst-level
B_file={}; B_state={}; B_burst=[]; B_center_rel_hz=[]; B_center_abs_hz=[];
B_cfo_coarse_hz=[]; B_cfo_fine_hz=[]; B_ngood_zc=[]; B_m1_b1=[]; B_m1_b2=[];
% 13-symbol block-level：NORMAL_64 / ALL_QPSK / MIXED
BL_file={}; BL_state={}; BL_burst=[]; BL_block=[]; BL_pattern={};
BL_nqpsk=[]; BL_n64=[]; BL_n16=[]; BL_m1rec=[];
% symbol-level
S_file={}; S_state={}; S_burst=[]; S_rel=[]; S_block=[]; S_m=[]; S_mod={};
S_evm=[]; S_kurt=[]; S_noiseVar=[]; S_cpCorr=[]; S_ncar=[];
% m1 detailed recurrence
M_file={}; M_state={}; M_burst=[]; M_block=[]; M_rel=[];
M_negI=[]; M_negQ=[]; M_posI=[]; M_posQ=[]; M_min=[];
% same-file repeated symbol consistency
R_file={}; R_state={}; R_rel=[]; R_m=[]; R_burstA=[]; R_burstB=[]; R_mismatch=[];
% global constellations
all64 = complex([]);
allM1 = complex([]);
%% ================================================================
% 4. File processing loop
% ================================================================
for fi=1:numel(matFiles)
    matFile = matFiles{fi};
    [~,baseName,~] = fileparts(matFile);
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s\n',fi,numel(matFiles),matFile);
    fprintf('============================================================\n');
    % ------------------------------------------------------------
    % 4.1 State labels for summaries only; not used for demodulation.
    % ------------------------------------------------------------
    if contains(matFile,char([hex2dec('98DE'),hex2dec('884C')]))
        stateLabel='flight';
    elseif contains(matFile,char([hex2dec('60AC'),hex2dec('505C')]))
        stateLabel='hover';
    elseif contains(matFile,char([hex2dec('9759'),hex2dec('6B62')]))
        stateLabel='static';
    else
        stateLabel='unknown';
    end
    % ------------------------------------------------------------
    % 4.2 Read IQ data, sampling frequency, and RF center.
    % ------------------------------------------------------------
    s=load(matFile);
    if isfield(s,'seg')
        x_raw=s.seg.iq;
        Fs_raw=double(s.seg.fs(1));
        Fc_raw=double(s.seg.fc(1));
    elseif isfield(s,'X')
        x_raw=s.X;
        if isfield(s,'Fs_out')
            Fs_raw=double(s.Fs_out(1));
        elseif isfield(s,'fs')
            Fs_raw=double(s.fs(1));
        else
            error('No sample rate in %s',matFile);
        end
        if isfield(s,'Fc_out')
            Fc_raw=double(s.Fc_out(1));
        else
            Fc_raw=0;
        end
    else
        error('Unsupported MAT IQ layout: %s',matFile);
    end
    if isreal(x_raw) && ismatrix(x_raw) && size(x_raw,2)==2
        x_raw=double(x_raw(:,1))+1j*double(x_raw(:,2));
    else
        x_raw=double(x_raw(:));
    end
    x_raw=x_raw(:);
    nRaw=numel(x_raw);
    durRaw=nRaw/Fs_raw;
    fprintf('  Fs/Fc       : %.6f MS/s / %.6f MHz\n',Fs_raw/1e6,Fc_raw/1e6);
    fprintf('  Samples     : %d (%.6f ms)\n',nRaw,1e3*durRaw);
    % The filename may say 200ms while the actual sample count represents only 20ms.
    % Lightweight spectrogram: FFT of approximately 300 short windows instead of a huge STFT array.
    if PLOT_OVERVIEW
        tfFFT=2048;
        if nRaw>=tfFFT
            fprintf('[PLOT] Figure 1: IQ spectrogram overview...\n');
            t0=tic;
            nTF=min(320,max(1,floor(nRaw/tfFFT)));
            tfBeg=round(linspace(1,nRaw-tfFFT+1,nTF));
            tfWin=hann(tfFFT,'periodic');
            tfDB=zeros(tfFFT,nTF,'single');
            for tfj=1:nTF
                tfIdx=tfBeg(tfj)+(0:tfFFT-1);
                tfVec=double(x_raw(tfIdx)).*tfWin;
                tfSpec=fftshift(fft(tfVec));
                tfDB(:,tfj)=single(20*log10(abs(tfSpec)+eps));
            end
            figure(1); clf;
            set(gcf,'Color','w','Name','Mini5Pro - Time-Frequency Preview','NumberTitle','off');
            tfFreq=(-tfFFT/2:tfFFT/2-1).'*Fs_raw/tfFFT;
            if Fc_raw~=0
                tfFreq=(tfFreq+Fc_raw)/1e6;
                tfLabel='RF frequency (MHz)';
            else
                tfFreq=tfFreq/1e6;
                tfLabel='Baseband offset (MHz)';
            end
            imagesc((tfBeg-1)/Fs_raw*1e3,tfFreq,tfDB);
            axis xy; xlabel('Time (ms)'); ylabel(tfLabel);
            title(sprintf('Mini5Pro IQ overview - %s',baseName),'Interpreter','none');
            colormap turbo; colorbar;
            drawnow;
            fprintf('[PLOT] Figure 1 ready (%.2f s).\n',toc(t0));
            clear tfDB tfVec tfSpec;
        end
    end
    % Initialize the constellation window before decoding; provide diagnostics if no bursts are found.
    if PLOT_CONSTELLATION
        figure(2); clf;
        set(gcf,'Color','w','Name','Mini5Pro - Live Constellations','NumberTitle','off', ...
            'Position',[100 100 1280 590]);
        axM1=subplot(1,2,1);
        axis(axM1,'equal'); grid(axM1,'on');
        xlim(axM1,[-1.4 1.4]); ylim(axM1,[-1.4 1.4]);
        xlabel(axM1,'I'); ylabel(axM1,'Q');
        title(axM1,'m1 QPSK - waiting for burst');
        axQAM=subplot(1,2,2);
        axis(axQAM,'equal'); grid(axQAM,'on');
        xlim(axQAM,[-1.4 1.4]); ylim(axQAM,[-1.4 1.4]);
        xlabel(axQAM,'I'); ylabel(axQAM,'Q');
        title(axQAM,'64QAM - waiting for burst');
        drawnow;
    end
    fileOut=fullfile(outputRoot,stateLabel,baseName);
    if ~exist(fileOut,'dir'), mkdir(fileOut); end
    % ------------------------------------------------------------
    % 4.3 Automatically find the 1201-tone video band center at the native FFT resolution.
    % ------------------------------------------------------------
    Nnative=round(NFFT*Fs_raw/FS_OUT);
    dec=round(Fs_raw/FS_OUT);
    if abs(Fs_raw/dec-FS_OUT)>1
        error('Fs %.9f MHz is not an integer multiple of 30.72 MHz',Fs_raw/1e6);
    end
    fprintf('[1/5] Searching occupied band...\n');
    phaseClock=tic;
    nf=floor(nRaw/Nnative);
    if nf<3
        warning('File too short for center search: %s',baseName);
        continue;
    end
    cand=zeros(nf,1);
    con=zeros(nf,1);
    cc=(700:Nnative-700).';
    Ps=zeros(Nnative,nf,'single');
    for ff=1:nf
        idx=(ff-1)*Nnative+(1:Nnative);
        Pw=abs(fftshift(fft(x_raw(idx),Nnative))).^2;
        Ps(:,ff)=single(Pw);
        cs=[0;cumsum(Pw)];
        mIn=(cs(cc+601)-cs(cc-600))/1201;
        mLo=(cs(cc-610)-cs(cc-640))/30;
        mHi=(cs(cc+641)-cs(cc+611))/30;
        ctr=10*log10(mIn./max(max(mLo,mHi),realmin));
        [con(ff),jj]=max(ctr);
        cand(ff)=cc(jj);
    end
    goodCenter=(con>10);
    if sum(goodCenter)<3
        tmp=sort(con);
        thr=tmp(max(1,ceil(0.70*numel(tmp))));
        goodCenter=(con>thr);
    end
    if ~any(goodCenter), goodCenter(:)=true; end
    cmed=median(cand(goodCenter));
    sel=goodCenter & abs(cand-cmed)<30;
    if ~any(sel), sel=goodCenter; end
    Psum=sum(double(Ps(:,sel)),2);
    clear Ps;
    cs=[0;cumsum(Psum)];
    mIn=(cs(cc+601)-cs(cc-600))/1201;
    mLo=(cs(cc-610)-cs(cc-640))/30;
    mHi=(cs(cc+641)-cs(cc+611))/30;
    ctr=10*log10(mIn./max(max(mLo,mHi),realmin));
    [boxCon,jj]=max(ctr);
    cb=cc(jj);
    fCenterHz=(cb-(Nnative/2+1))*Fs_raw/Nnative;
    fprintf('  Band center : %+.6f MHz rel Fc | box contrast %.2f dB (%.2f s)\n',fCenterHz/1e6,boxCon,toc(phaseClock));
    % ------------------------------------------------------------
    % 4.4 DDC + anti-alias LPF + decimate -> 30.72 MS/s
    % ------------------------------------------------------------
    fprintf('[2/5] DDC and 601-tap FIR: FFT convolution...\n');
    phaseClock=tic;
    nraw=(0:nRaw-1).';
    xm=x_raw.*exp(-1j*2*pi*fCenterHz/Fs_raw*nraw);
    clear x_raw nraw;
    bLP=fir1(600,9.45e6/(Fs_raw/2),kaiser(601,9));
    % fftfilt and filter(bLP,1,xm) implement the same causal FIR convolution.
    % Some older MATLAB versions require double input for fftfilt; convert explicitly.
    yf=fftfilt(double(bLP(:)),double(xm(:)));
    clear xm;
    % Compensate the FIR group delay of 300 samples, retaining the original time origin.
    yf=[yf(301:end);complex(zeros(300,1))];
    y0=yf(1:dec:end);
    clear yf;
    fprintf('[2/5] DDC/FIR finished in %.2f s; output samples=%d\n',toc(phaseClock),numel(y0));
    y0=y0-mean(y0);
    y0=y0/sqrt(mean(abs(y0).^2)+eps);
    ny=(0:numel(y0)-1).';
    % ------------------------------------------------------------
    % 4.5 Coarse fractional CFO estimation from high CP correlation.
    % ------------------------------------------------------------
    fprintf('[3/5] Fractional CFO from CP...\n');
    phaseClock=tic;
    c=y0(1:end-NFFT).*conj(y0(1+NFFT:end));
    % Replace three 128-point convolutions with prefix sums; preserve valid-convolution indexing.
    cpSum=[0;cumsum(c)];
    cs128=cpSum(129:end)-cpSum(1:end-128);
    cpEnergy1=abs(y0(1:end-NFFT)).^2;
    cpEnergy2=abs(y0(1+NFFT:end)).^2;
    cpCum1=[0;cumsum(cpEnergy1)];
    cpCum2=[0;cumsum(cpEnergy2)];
    e1=cpCum1(129:end)-cpCum1(1:end-128);
    e2=cpCum2(129:end)-cpCum2(1:end-128);
    clear c cpSum cpEnergy1 cpEnergy2 cpCum1 cpCum2;
    e1=max(e1,0); e2=max(e2,0);
    rr=abs(cs128)./sqrt(e1.*e2+eps);
    if any(rr>0.95)
        fFrac=-angle(sum(cs128(rr>0.95)))*FS_OUT/(2*pi*NFFT);
    else
        fFrac=0;
    end
    y0=y0.*exp(-1j*2*pi*fFrac/FS_OUT*ny);
    fprintf('  Fraction CFO: %+.2f Hz (%.2f s)\n',fFrac,toc(phaseClock));
    % ------------------------------------------------------------
    % 4.6 Joint integer CFO (+/-8 bins) selection using ZC and active-band edge contrast.
    % ------------------------------------------------------------
    fprintf('[4/5] Integer CFO sweep -8..+8 bins...\n');
    phaseClock=tic;
    bestScore=-inf;
    kcBest=0;
    Ny0=numel(y0);
    nValid0=Ny0-NFFT+1;
    % The ZC matched-filter FFT does not need to be recomputed for every candidate.
    hMatchFFT=conj(fft(tZ,Ny0));
    % Digital frequency shifting changes phase only, so sliding signal energy is shared among all 17 candidates.
    energyCum=[0;cumsum(abs(y0).^2)];
    pZ0=max(energyCum(NFFT+1:end)-energyCum(1:end-NFFT),0);
    clear energyCum;
    for kc=-8:8
        yk=y0.*exp(-1j*2*pi*(kc*SCS)/FS_OUT*ny);
        Ny=Ny0;
        ccZ=ifft(fft(yk).*hMatchFFT);
        rhoK=abs(ccZ(1:nValid0))./sqrt(pZ0+eps);
        % NMS peak picker：threshold 0.8, +/-1500
        tmpR=rhoK;
        pkK=[];
        while numel(pkK)<400
            [vv,iiPeak]=max(tmpR);
            if vv<INTEGER_ZC_RHO_MIN, break; end
            pkK(end+1,1)=iiPeak; %#ok<SAGROW>
            a=max(1,iiPeak-1500);
            b=min(numel(tmpR),iiPeak+1500);
            tmpR(a:b)=0;
        end
        pkK=sort(pkK);
        if numel(pkK)<3
            fprintf('  bin %+d / 8: ZC peaks=%d, skipped\n',kc,numel(pkK));
            drawnow limitrate;
            continue;
        end
        Pw=zeros(NFFT,1);
        for ip=1:numel(pkK)
            pp=pkK(ip);
            if pp+NFFT-1<=Ny
                Pw=Pw+abs(fft(yk(pp:pp+NFFT-1),NFFT)).^2;
            end
        end
        idxInL=mod((-600:-595).',NFFT)+1;
        idxInR=mod((595:600).',NFFT)+1;
        idxOutL=mod((-606:-601).',NFFT)+1;
        idxOutR=mod((601:606).',NFFT)+1;
        inE=min(mean(Pw(idxInL)),mean(Pw(idxInR)));
        outE=max(mean(Pw(idxOutL)),mean(Pw(idxOutR)));
        sc=10*log10(inE/(outE+eps));
        fprintf('  bin %+d / 8: ZC peaks=%d, edge contrast=%.2f dB\n',kc,numel(pkK),sc);
        drawnow limitrate;
        if sc>bestScore
            bestScore=sc;
            kcBest=kc;
        end
    end
    fprintf('[4/5] Integer CFO sweep finished in %.2f s\n',toc(phaseClock));
    clear yk rhoK tmpR ccZ pkK Pw;
    if ~isfinite(bestScore)
        warning('No integer CFO candidate reached ZC threshold; fallback kc=0.');
    end
    if bestScore<EDGE_CONTRAST_MINDB
        warning('Integer-bin decision weak in %s: edge contrast %.2f dB',baseName,bestScore);
    end
    y=y0.*exp(-1j*2*pi*(kcBest*SCS)/FS_OUT*ny);
    clear y0;
    centerRelHz=fCenterHz+fFrac+kcBest*SCS;
    centerAbsHz=Fc_raw+centerRelHz;
    fprintf('  Integer CFO : %+d bins = %+.3f kHz | edge %.2f dB\n',kcBest,kcBest*SCS/1e3,bestScore);
    fprintf('  Final center: %.6f MHz absolute (%+.6f MHz rel Fc)\n',centerAbsHz/1e6,centerRelHz/1e6);
    % ------------------------------------------------------------
    % 4.7 ZC peaks and complete burst structure identification.
    % ------------------------------------------------------------
    fprintf('[5/5] Final ZC correlation and burst decode...\n');
    phaseClock=tic;
    Ny=numel(y);
    % Reuse hMatchFFT / pZ0: the final signal differs from y0 only by an integer-bin frequency rotation.
    ccZ=ifft(fft(y).*hMatchFFT);
    rho=abs(ccZ(1:nValid0))./sqrt(pZ0+eps);
    clear ccZ hMatchFFT pZ0;
    tmpR=rho;
    pk=[];
    while numel(pk)<400
        [vv,iiPeak]=max(tmpR);
        if vv<ZC_RHO_MIN, break; end
        pk(end+1,1)=iiPeak; %#ok<SAGROW>
        a=max(1,iiPeak-1500);
        b=min(numel(tmpR),iiPeak+1500);
        tmpR(a:b)=0;
    end
    pk=sort(pk);
    bursts={};
    for ii=1:max(0,numel(pk)-4)
        ddp=diff(pk(ii:ii+4));
        if ~all(abs(ddp(:).'-[6848 15360 13168 15360])<=3)
            continue;
        end
        z1=pk(ii+1);
        bestS=-inf;
        b3=z1;
        for candB3=(z1-3):(z1+3)
            testIdx=candB3+[0 15360 28528 43888];
            if min(testIdx)<1 || max(testIdx)>numel(rho), continue; end
            ss=sum(rho(testIdx));
            if ss>bestS
                bestS=ss;
                b3=candB3;
            end
        end
        zz=[b3-6848 b3 b3+15360 b3+28528 b3+43888];
        B=zeros(27,1);
        B(1)=zz(1);
        for rel=1:26
            if rel<=9
                B(rel+1)=b3+2192*(rel-3)+16*(rel>=8);
            elseif rel<=15
                B(rel+1)=zz(3)+2192*(rel-10)+16*(rel>=15);
            elseif rel<=22
                B(rel+1)=zz(4)+2192*(rel-16)+16*(rel>=22);
            else
                B(rel+1)=zz(5)+2192*(rel-23);
            end
        end
        if zz(1)>300 && B(27)+NFFT-1<=Ny
            isDup=false;
            for jj=1:numel(bursts)
                if abs(bursts{jj}(1)-zz(1))<100
                    isDup=true;
                    break;
                end
            end
            if ~isDup
                bursts{end+1,1}=zz; %#ok<SAGROW>
            end
        end
    end
    fprintf('  ZC peaks     : %d | complete bursts %d (%.2f s)\n',numel(pk),numel(bursts),toc(phaseClock));
    % Always include file-level summaries even when no complete burst is detected.
    F_file{end+1,1}=baseName; %#ok<SAGROW>
    F_state{end+1,1}=stateLabel; %#ok<SAGROW>
    F_fs(end+1,1)=Fs_raw; %#ok<SAGROW>
    F_fc(end+1,1)=Fc_raw; %#ok<SAGROW>
    F_nsamp(end+1,1)=nRaw; %#ok<SAGROW>
    F_duration(end+1,1)=durRaw; %#ok<SAGROW>
    F_center_rel_hz(end+1,1)=centerRelHz; %#ok<SAGROW>
    F_center_abs_hz(end+1,1)=centerAbsHz; %#ok<SAGROW>
    F_edge_db(end+1,1)=bestScore; %#ok<SAGROW>
    F_nburst(end+1,1)=numel(bursts); %#ok<SAGROW>
    if isempty(bursts)
        warning('No complete Mini5Pro burst found in %s; no constellation can be reconstructed.',baseName);
        if PLOT_CONSTELLATION
            title(axM1,'No complete ZC burst - no QPSK');
            title(axQAM,'No complete ZC burst - no QAM');
            drawnow;
        end
        continue;
    end
    % ------------------------------------------------------------
    % 4.8 Current-file recovery buffers
    % ------------------------------------------------------------
    nB=numel(bursts);
    fileZ=complex(nan(nB,27,L,'single'));
    fileH=complex(nan(nB,27,L,'single'));
    fileLev=zeros(nB,27,L,2,'int8');
    fileMod=zeros(nB,27,'uint8');      % 0=ZC/none, 1=QPSK, 2=16QAM, 3=64QAM
    fileEvm=nan(nB,27,'single');
    fileCp=nan(nB,27,'single');
    file64=complex([]);
    fileM1=complex([]);
    % ============================================================
    % 5. burst loop
    % ============================================================
    for bi=1:nB
        fprintf('[BURST %d/%d] start\n', bi, nB);
        burstClock = tic;
        currentBurstCell = bursts(bi);
        zz = currentBurstCell{1};
        b3=zz(2);
        B=zeros(27,1);
        B(1)=zz(1);
        for rel=1:26
            if rel<=9
                B(rel+1)=b3+2192*(rel-3)+16*(rel>=8);
            elseif rel<=15
                B(rel+1)=zz(3)+2192*(rel-10)+16*(rel>=15);
            elseif rel<=22
                B(rel+1)=zz(4)+2192*(rel-16)+16*(rel>=22);
            else
                B(rel+1)=zz(5)+2192*(rel-23);
            end
        end
        lo=max(1,B(1)-3000);
        hi=min(Ny,B(27)+NFFT+100);
        seg=y(lo:hi);
        nn=(lo-1:hi-1).';
        % --------------------------------------------------------
        % 5.1 Coarse CFO from the last 128 CP samples across 27 symbol positions.
        % --------------------------------------------------------
        acc=0;
        for rel=0:26
            p0=B(rel+1)-lo+1;
            a1=p0-128;
            a2=p0-1;
            if a1<1 || a2+NFFT>numel(seg), continue; end
            acc=acc+sum(seg(a1:a2).*conj(seg(a1+NFFT:a2+NFFT)));
        end
        cfoCoarse=-angle(acc)*FS_OUT/(2*pi*NFFT);
        seg=seg.*exp(-1j*2*pi*cfoCoarse/FS_OUT*nn);
        % --------------------------------------------------------
        % 5.2 Fine CFO from common phase slope of four ZC symbols; reject incoherent anchors.
        % --------------------------------------------------------
        H0=cell(4,1);
        for aa=1:4
            rel=CFO_ZC_REL(aa);
            p=B(rel+1)-lo+1;
            Y=fft(seg(p:p+NFFT-1),NFFT);
            Y=Y(KIDX);
            H0{aa}=Y.*conj(ZC);
        end
        Cmat=zeros(4,4);
        for aa=1:4
            for bb=1:4
                if aa==bb, continue; end
                Cmat(aa,bb)=abs(sum(conj(H0{aa}).*H0{bb})) / ...
                    (sqrt(sum(abs(H0{aa}).^2))*sqrt(sum(abs(H0{bb}).^2))+eps);
            end
        end
        goodIdx=find(max(Cmat,[],2)>0.98).';
        if numel(goodIdx)<2, goodIdx=1:4; end
        goodRel=CFO_ZC_REL(goodIdx);
        th=zeros(1,numel(goodRel));
        tt=zeros(1,numel(goodRel));
        Href=H0{goodIdx(1)};
        for aa=1:numel(goodRel)
            th(aa)=angle(sum(conj(Href).*H0{goodIdx(aa)}));
            tt(aa)=B(goodRel(aa)+1);
        end
        pfc=polyfit(tt,unwrap(th),1);
        cfoFine=pfc(1)*FS_OUT/(2*pi);
        seg=seg.*exp(-1j*pfc(1)*nn);
        % --------------------------------------------------------
        % 5.3 ZC least-squares channel and common-phase trend compensation.
        % --------------------------------------------------------
        Hr=cell(27,1);
        for aa=1:numel(ANCHOR_REL)
            rel=ANCHOR_REL(aa);
            p=B(rel+1)-lo+1;
            Y=fft(seg(p:p+NFFT-1),NFFT);
            Y=Y(KIDX);
            Hr{rel+1}=Y.*conj(ZC);
        end
        th2=zeros(1,numel(goodRel));
        Href=Hr{goodRel(1)+1};
        for aa=1:numel(goodRel)
            th2(aa)=angle(sum(conj(Href).*Hr{goodRel(aa)+1}));
        end
        pf=polyfit(tt,unwrap(th2),1);
        Hd=cell(27,1);
        for aa=1:numel(ANCHOR_REL)
            rel=ANCHOR_REL(aa);
            Hd{rel+1}=Hr{rel+1}.*exp(-1j*polyval(pf,B(rel+1)));
        end
        % Current burst m1 self-validation buffer.
        chkBlock=nan(1,2);
        % --------------------------------------------------------
        % 5.4 Demodulate data positions rel=1..26.
        % --------------------------------------------------------
        for rel=1:26
            if any(rel==ANCHOR_REL), continue; end
            left=goodRel(goodRel<rel);
            right=goodRel(goodRel>rel);
            if ~isempty(left) && ~isempty(right)
                aL=max(left); aR=min(right);
            elseif ~isempty(right)
                aL=right(1);
                aR=right(min(2,numel(right)));
            elseif ~isempty(left)
                aL=max(left); aR=aL;
            else
                continue;
            end
            if aL==aR
                H=Hd{aL+1};
            else
                al=(B(rel+1)-B(aL+1))/(B(aR+1)-B(aL+1));
                H=(1-al)*Hd{aL+1}+al*Hd{aR+1};
            end
            H=H.*exp(1j*polyval(pf,B(rel+1)));
            p=B(rel+1)-lo+1;
            Y=fft(seg(p:p+NFFT-1),NFFT);
            Y=Y(KIDX);
            Z=Y./(H+eps);
            zv=Z(NODC);
            ku=mean(abs(zv).^4)/(mean(abs(zv).^2)^2+eps);
            % ---------- modulation / EVM ----------
            z0=zv/sqrt(mean(abs(zv).^2)+eps);
            % rel=1 and rel=14 are known m1 QPSK positions; do not let noise-corrupted kurtosis
            % misclassify them as QAM and leave the original m1 constellation empty.
            if rel==1 || rel==14 || ku<1.15
                modu='QPSK';
                Cuse=CQ;
                modID=1;
            else
                % 64QAM EVM candidate
                a64=1;
                for it=1:6
                    DD=abs(z0/a64-C64.').^2;
                    [~,id]=min(DD,[],2);
                    q=C64(id);
                    a64=real(sum(conj(q).*z0))/(sum(abs(q).^2)+eps);
                end
                zn64=z0/a64;
                DD=abs(zn64-C64.').^2;
                [~,id]=min(DD,[],2);
                q=C64(id);
                ev64=sqrt(mean(abs(zn64-q).^2)/(mean(abs(q).^2)+eps));
                % 16QAM EVM candidate
                a16=1;
                for it=1:6
                    DD=abs(z0/a16-C16.').^2;
                    [~,id]=min(DD,[],2);
                    q=C16(id);
                    a16=real(sum(conj(q).*z0))/(sum(abs(q).^2)+eps);
                end
                zn16=z0/a16;
                DD=abs(zn16-C16.').^2;
                [~,id]=min(DD,[],2);
                q=C16(id);
                ev16=sqrt(mean(abs(zn16-q).^2)/(mean(abs(q).^2)+eps));
                if ev16<ev64
                    modu='16QAM';
                    Cuse=C16;
                    modID=2;
                else
                    modu='64QAM';
                    Cuse=C64;
                    modID=3;
                end
            end
            aScale=1;
            for it=1:6
                DD=abs(z0/aScale-Cuse.').^2;
                [~,id]=min(DD,[],2);
                q=Cuse(id);
                aScale=real(sum(conj(q).*z0))/(sum(abs(q).^2)+eps);
            end
            znNoDC=z0/aScale;
            DD=abs(znNoDC-Cuse.').^2;
            [~,id]=min(DD,[],2);
            q=Cuse(id);
            ev=sqrt(mean(abs(znNoDC-q).^2)/(mean(abs(q).^2)+eps));
            % Apply the same normalization to all 1201 bins, including DC.
            zn=Z/(sqrt(mean(abs(zv).^2)+eps)*aScale);
            nv=ev^2;
            fileZ(bi,rel+1,:)=reshape(single(zn),1,1,L);
            fileH(bi,rel+1,:)=reshape(single(H),1,1,L);
            fileMod(bi,rel+1)=uint8(modID);
            fileEvm(bi,rel+1)=single(ev);
            % geometry level
            if modID==3
                Il=2*round((real(zn)*sqrt(42)-1)/2)+1;
                Ql=2*round((imag(zn)*sqrt(42)-1)/2)+1;
                Il=min(max(Il,-7),7);
                Ql=min(max(Ql,-7),7);
                fileLev(bi,rel+1,:,1)=reshape(int8(Il),1,1,L);
                fileLev(bi,rel+1,:,2)=reshape(int8(Ql),1,1,L);
                file64=[file64;zn(NODC)]; %#ok<AGROW>
                all64=[all64;zn(NODC)]; %#ok<AGROW>
            elseif modID==1
                Il=sign(real(zn)); Ql=sign(imag(zn));
                Il(Il==0)=1; Ql(Ql==0)=1;
                fileLev(bi,rel+1,:,1)=reshape(int8(Il),1,1,L);
                fileLev(bi,rel+1,:,2)=reshape(int8(Ql),1,1,L);
                if rel==1 || rel==14
                    fileM1=[fileM1;zn(NODC)]; %#ok<AGROW>
                    allM1=[allM1;zn(NODC)]; %#ok<AGROW>
                end
            end
            % CP correlation
            Lcp=CP_N;
            if any(rel==LONG_CP_REL), Lcp=CP_L; end
            p0=B(rel+1)-lo+1;
            cpc=NaN;
            if p0-Lcp>=1 && p0-1+NFFT<=numel(seg)
                aa=seg(p0-Lcp:p0-1);
                bb=seg(p0-Lcp+NFFT:p0-1+NFFT);
                cpc=abs(sum(conj(bb).*aa))/(sqrt(sum(abs(aa).^2))*sqrt(sum(abs(bb).^2))+eps);
            end
            fileCp(bi,rel+1)=single(cpc);
            blockNo=1+(rel>=14);
            mNo=mod(rel-1,13)+1;
            S_file{end+1,1}=baseName; %#ok<SAGROW>
            S_state{end+1,1}=stateLabel; %#ok<SAGROW>
            S_burst(end+1,1)=bi-1; %#ok<SAGROW>
            S_rel(end+1,1)=rel; %#ok<SAGROW>
            S_block(end+1,1)=blockNo; %#ok<SAGROW>
            S_m(end+1,1)=mNo; %#ok<SAGROW>
            S_mod{end+1,1}=modu; %#ok<SAGROW>
            S_evm(end+1,1)=ev; %#ok<SAGROW>
            S_kurt(end+1,1)=ku; %#ok<SAGROW>
            S_noiseVar(end+1,1)=nv; %#ok<SAGROW>
            S_cpCorr(end+1,1)=cpc; %#ok<SAGROW>
            S_ncar(end+1,1)=1200; %#ok<SAGROW>
        end
        % --------------------------------------------------------
        % 5.5 m1 recurrence check: rel1=P131, rel14=P95.
        % --------------------------------------------------------
        for jj=1:2
            if jj==1
                rel=1; taps=TAPS131;
            else
                rel=14; taps=TAPS95;
            end
            ztmp=squeeze(fileZ(bi,rel+1,:));
            if all(isnan(ztmp)), continue; end
            acc4=nan(1,4);
            halfSign=[-1 1];
            ccnt=0;
            for hh=1:2
                selHalf=(sign(K)==halfSign(hh));
                zr=ztmp(selHalf);
                for iq=1:2
                    ccnt=ccnt+1;
                    if iq==1
                        bits=double(real(zr)<0);
                    else
                        bits=double(imag(zr)<0);
                    end
                    Lm=max(taps);
                    nb=numel(bits);
                    if nb<=Lm
                        acc4(ccnt)=NaN;
                    else
                        vv=bits(Lm+1:nb);
                        for tp=taps
                            vv=mod(vv+bits(Lm+1-tp:nb-tp),2);
                        end
                        acc4(ccnt)=1-mean(vv);
                    end
                end
            end
            rmin=min(acc4);
            chkBlock(jj)=rmin;
            M_file{end+1,1}=baseName; %#ok<SAGROW>
            M_state{end+1,1}=stateLabel; %#ok<SAGROW>
            M_burst(end+1,1)=bi-1; %#ok<SAGROW>
            M_block(end+1,1)=jj; %#ok<SAGROW>
            M_rel(end+1,1)=rel; %#ok<SAGROW>
            M_negI(end+1,1)=acc4(1); %#ok<SAGROW>
            M_negQ(end+1,1)=acc4(2); %#ok<SAGROW>
            M_posI(end+1,1)=acc4(3); %#ok<SAGROW>
            M_posQ(end+1,1)=acc4(4); %#ok<SAGROW>
            M_min(end+1,1)=rmin; %#ok<SAGROW>
        end
        % --------------------------------------------------------
        % 5.6 Modulation-pattern classification for two 13-symbol blocks.
        % Each block has 13 rel positions; two carry ZC, leaving 11 data symbols.
        % NORMAL_64: m1 is QPSK and the other 10 data symbols are 64QAM.
        % ALL_QPSK: all 11 data symbols are QPSK.
        % --------------------------------------------------------
        for blockNo=1:2
            if blockNo==1
                rrBlock=1:13;
            else
                rrBlock=14:26;
            end
            rrData=rrBlock(~ismember(rrBlock,ANCHOR_REL));
            mids=double(fileMod(bi,rrData+1));
            nq=sum(mids==1); n16=sum(mids==2); n64=sum(mids==3);
            if nq==11 && n16==0 && n64==0
                patt='ALL_QPSK';
            elseif nq==1 && n16==0 && n64==10
                patt='NORMAL_64';
            else
                patt='MIXED';
            end
            BL_file{end+1,1}=baseName; %#ok<SAGROW>
            BL_state{end+1,1}=stateLabel; %#ok<SAGROW>
            BL_burst(end+1,1)=bi-1; %#ok<SAGROW>
            BL_block(end+1,1)=blockNo; %#ok<SAGROW>
            BL_pattern{end+1,1}=patt; %#ok<SAGROW>
            BL_nqpsk(end+1,1)=nq; %#ok<SAGROW>
            BL_n64(end+1,1)=n64; %#ok<SAGROW>
            BL_n16(end+1,1)=n16; %#ok<SAGROW>
            BL_m1rec(end+1,1)=chkBlock(blockNo); %#ok<SAGROW>
        end
        B_file{end+1,1}=baseName; %#ok<SAGROW>
        B_state{end+1,1}=stateLabel; %#ok<SAGROW>
        B_burst(end+1,1)=bi-1; %#ok<SAGROW>
        B_center_rel_hz(end+1,1)=centerRelHz; %#ok<SAGROW>
        B_center_abs_hz(end+1,1)=centerAbsHz; %#ok<SAGROW>
        B_cfo_coarse_hz(end+1,1)=cfoCoarse; %#ok<SAGROW>
        B_cfo_fine_hz(end+1,1)=cfoFine; %#ok<SAGROW>
        B_ngood_zc(end+1,1)=numel(goodRel); %#ok<SAGROW>
        B_m1_b1(end+1,1)=chkBlock(1); %#ok<SAGROW>
        B_m1_b2(end+1,1)=chkBlock(2); %#ok<SAGROW>
        fprintf('  burst %02d: CFO %+8.1f %+7.1f Hz | goodZC=%s | m1 %.4f / %.4f\n', ...
            bi-1,cfoCoarse,cfoFine,mat2str(goodRel),chkBlock(1),chkBlock(2));
        if any(chkBlock<M1_WARN)
            warning('m1 recurrence below %.4f in %s burst %d',M1_WARN,baseName,bi-1);
        end
        % Update live constellations during burst processing; no need to await MAT/CSV export.
        if PLOT_CONSTELLATION && (bi==1 || mod(bi,CONSTELLATION_UPDATE_EVERY)==0 || bi==nB)
            cla(axM1); cla(axQAM);
            if ~isempty(fileM1)
                qDraw=fileM1;
                if numel(qDraw)>MAX_CONSTELLATION_POINTS
                    pick=round(linspace(1,numel(qDraw),MAX_CONSTELLATION_POINTS));
                    qDraw=qDraw(pick);
                end
                plot(axM1,real(qDraw),imag(qDraw),'.','MarkerSize',3);
            else
                text(axM1,0,0,'No valid m1 QPSK yet','HorizontalAlignment','center');
            end
            grid(axM1,'on'); axis(axM1,'equal');
            xlim(axM1,[-1.4 1.4]); ylim(axM1,[-1.4 1.4]);
            xlabel(axM1,'I'); ylabel(axM1,'Q');
            title(axM1,sprintf('m1 QPSK | burst %d/%d | n=%d',bi,nB,numel(fileM1)));
            if ~isempty(file64)
                zDraw=file64;
                if numel(zDraw)>MAX_CONSTELLATION_POINTS
                    pick=round(linspace(1,numel(zDraw),MAX_CONSTELLATION_POINTS));
                    zDraw=zDraw(pick);
                end
                plot(axQAM,real(zDraw),imag(zDraw),'.','MarkerSize',3);
                hold(axQAM,'on');
                plot(axQAM,real(C64),imag(C64),'kx','MarkerSize',5,'LineWidth',1);
                hold(axQAM,'off');
            else
                text(axQAM,0,0,'No 64QAM decoded yet (possible ALL_QPSK)', ...
                    'HorizontalAlignment','center');
            end
            grid(axQAM,'on'); axis(axQAM,'equal');
            xlim(axQAM,[-1.4 1.4]); ylim(axQAM,[-1.4 1.4]);
            xlabel(axQAM,'I'); ylabel(axQAM,'Q');
            title(axQAM,sprintf('64QAM | burst %d/%d | n=%d',bi,nB,numel(file64)));
            drawnow limitrate;
        end
        fprintf('[BURST %d/%d] Finished (%.2f s), QPSK=%d, 64QAM=%d\n', ...
            bi,nB,toc(burstClock),numel(fileM1),numel(file64));
    end
    %% ============================================================
    % 6. Same-position 64QAM repeatability analysis within the current file.
    % ============================================================
    for rel=1:26
        ids=find(fileMod(:,rel+1)==3);
        if numel(ids)<2, continue; end
        for ia=1:numel(ids)-1
            for ib=ia+1:numel(ids)
                ba=ids(ia); bb=ids(ib);
                A=squeeze(fileLev(ba,rel+1,NODC,:));
                C=squeeze(fileLev(bb,rel+1,NODC,:));
                mismatch=mean(any(A~=C,2));
                R_file{end+1,1}=baseName; %#ok<SAGROW>
                R_state{end+1,1}=stateLabel; %#ok<SAGROW>
                R_rel(end+1,1)=rel; %#ok<SAGROW>
                R_m(end+1,1)=mod(rel-1,13)+1; %#ok<SAGROW>
                R_burstA(end+1,1)=ba-1; %#ok<SAGROW>
                R_burstB(end+1,1)=bb-1; %#ok<SAGROW>
                R_mismatch(end+1,1)=mismatch; %#ok<SAGROW>
            end
        end
    end
    %% ============================================================
    % 7. Per-file outputs
    % ============================================================
    if SAVE_PER_FILE_MAT
        Mini5Pro.file=baseName;
        Mini5Pro.state=stateLabel;
        Mini5Pro.fs_raw=Fs_raw;
        Mini5Pro.fc_raw=Fc_raw;
        Mini5Pro.center_rel_hz=centerRelHz;
        Mini5Pro.center_abs_hz=centerAbsHz;
        Mini5Pro.K=K;
        Mini5Pro.Z=fileZ;
        Mini5Pro.H=fileH;
        Mini5Pro.level=fileLev;
        Mini5Pro.modulation=fileMod;
        Mini5Pro.evm=fileEvm;
        Mini5Pro.cpCorr=fileCp;
        Mini5Pro.modulation_legend='0=ZC/none,1=QPSK,2=16QAM,3=64QAM';
        Mini5Pro.note='64QAM I/Q level is geometry only; final DJI bit labeling not proven.';
        save(fullfile(fileOut,[baseName '_Mini5Pro_V4_recovered.mat']),'Mini5Pro','-v7.3');
        clear Mini5Pro;
    end
    if SAVE_HARD_LEVEL_CSV
        hardRows=[];
        for bi=1:nB
            for rel=1:26
                if fileMod(bi,rel+1)~=3, continue; end
                Il=double(squeeze(fileLev(bi,rel+1,NODC,1)));
                Ql=double(squeeze(fileLev(bi,rel+1,NODC,2)));
                knd=double(K(NODC));
                nr=numel(knd);
                hardRows=[hardRows; ...
                    [repmat(bi-1,nr,1),repmat(rel,nr,1),knd,Il,Ql]]; %#ok<AGROW>
            end
        end
        if ~isempty(hardRows)
            Thard=array2table(hardRows,'VariableNames',{'burst','rel','K','I_level','Q_level'});
            writetable(Thard,fullfile(fileOut,[baseName '_64QAM_hard_levels.csv']));
        end
    end
    if PLOT_CONSTELLATION && ~isempty(file64)
        figure(3); clf; set(gcf,'Color','w','Name','Mini5Pro - File 64QAM');
        qPlot=file64; if numel(qPlot)>MAX_CONSTELLATION_POINTS, qPlot=qPlot(round(linspace(1,numel(qPlot),MAX_CONSTELLATION_POINTS))); end
        scatter(real(qPlot),imag(qPlot),2,'.');
        axis equal; grid on;
        xlim([-1.4 1.4]); ylim([-1.4 1.4]);
        xlabel('I'); ylabel('Q');
        title(sprintf('Mini5Pro V4 physical-layer 64QAM - %s',baseName),'Interpreter','none');
        saveas(gcf,fullfile(fileOut,[baseName '_64QAM_constellation.png']));
        drawnow;  % Keep the constellation window open for inspection (the old version closed it).
    end
    if PLOT_CONSTELLATION && ~isempty(fileM1)
        figure(4); clf; set(gcf,'Color','w','Name','Mini5Pro - File m1 QPSK');
        qPlot=fileM1; if numel(qPlot)>MAX_CONSTELLATION_POINTS, qPlot=qPlot(round(linspace(1,numel(qPlot),MAX_CONSTELLATION_POINTS))); end
        scatter(real(qPlot),imag(qPlot),3,'.');
        axis equal; grid on;
        xlim([-1.2 1.2]); ylim([-1.2 1.2]);
        xlabel('I'); ylabel('Q');
        title(sprintf('Mini5Pro V4 m1 QPSK - %s',baseName),'Interpreter','none');
        saveas(gcf,fullfile(fileOut,[baseName '_m1_QPSK_constellation.png']));
        drawnow;  % Keep the constellation window open for inspection (the old version closed it).
    end
end
%% ================================================================
% 8. Global CSV summaries
% ================================================================
if ~isempty(F_file)
    Tfile=table(F_file,F_state,F_fs,F_fc,F_nsamp,F_duration,F_center_rel_hz,F_center_abs_hz,F_edge_db,F_nburst, ...
        'VariableNames',{'file','state','fs_raw_hz','fc_raw_hz','samples','duration_s', ...
        'center_rel_hz','center_abs_hz','edge_contrast_db','n_burst'});
    writetable(Tfile,fullfile(outputRoot,'Mini5Pro_V4_file_summary.csv'));
end
if ~isempty(B_file)
    Tburst=table(B_file,B_state,B_burst,B_center_rel_hz,B_center_abs_hz, ...
        B_cfo_coarse_hz,B_cfo_fine_hz,B_ngood_zc,B_m1_b1,B_m1_b2, ...
        'VariableNames',{'file','state','burst','center_rel_hz','center_abs_hz', ...
        'cfo_coarse_hz','cfo_fine_hz','n_good_zc','m1_block1','m1_block2'});
    writetable(Tburst,fullfile(outputRoot,'Mini5Pro_V4_burst_summary.csv'));
end
if ~isempty(BL_file)
    Tblock=table(BL_file,BL_state,BL_burst,BL_block,BL_pattern,BL_nqpsk,BL_n64,BL_n16,BL_m1rec, ...
        'VariableNames',{'file','state','burst','block','pattern','n_qpsk','n_64qam','n_16qam','m1_recurrence'});
    writetable(Tblock,fullfile(outputRoot,'Mini5Pro_V4_block_summary.csv'));
end
if ~isempty(S_file)
    Tsym=table(S_file,S_state,S_burst,S_rel,S_block,S_m,S_mod,S_evm,S_kurt,S_noiseVar,S_cpCorr,S_ncar, ...
        'VariableNames',{'file','state','burst','rel','block','m','modulation','evm', ...
        'kurtosis','noise_var','cp_corr','ncar'});
    writetable(Tsym,fullfile(outputRoot,'Mini5Pro_V4_symbol_summary.csv'));
end
if ~isempty(M_file)
    Tm1=table(M_file,M_state,M_burst,M_block,M_rel,M_negI,M_negQ,M_posI,M_posQ,M_min, ...
        'VariableNames',{'file','state','burst','block','rel','negI','negQ','posI','posQ','min_recurrence'});
    writetable(Tm1,fullfile(outputRoot,'Mini5Pro_V4_m1_recurrence.csv'));
end
if ~isempty(R_file)
    Trepeat=table(R_file,R_state,R_rel,R_m,R_burstA,R_burstB,R_mismatch, ...
        'VariableNames',{'file','state','rel','m','burst_a','burst_b','symbol_mismatch_rate'});
    writetable(Trepeat,fullfile(outputRoot,'Mini5Pro_V4_same_position_repeat.csv'));
    % Group by flight state and rel to identify stable/filler positions.
    RR_state={}; RR_rel=[]; RR_m=[]; RR_pairs=[]; RR_med=[]; RR_min=[]; RR_le1=[]; RR_ge90=[];
    rStates=unique(Trepeat.state,'stable');
    for is=1:numel(rStates)
        st=rStates{is};
        for rel=1:26
            sel=strcmp(Trepeat.state,st) & (Trepeat.rel==rel);
            if ~any(sel), continue; end
            vv=Trepeat.symbol_mismatch_rate(sel);
            RR_state{end+1,1}=st; %#ok<SAGROW>
            RR_rel(end+1,1)=rel; %#ok<SAGROW>
            RR_m(end+1,1)=mod(rel-1,13)+1; %#ok<SAGROW>
            RR_pairs(end+1,1)=numel(vv); %#ok<SAGROW>
            RR_med(end+1,1)=median(vv); %#ok<SAGROW>
            RR_min(end+1,1)=min(vv); %#ok<SAGROW>
            RR_le1(end+1,1)=sum(vv<=1/1200); %#ok<SAGROW>
            RR_ge90(end+1,1)=sum(vv>=0.90); %#ok<SAGROW>
        end
    end
    TrepeatRel=table(RR_state,RR_rel,RR_m,RR_pairs,RR_med,RR_min,RR_le1,RR_ge90, ...
        'VariableNames',{'state','rel','m','n_pairs','mismatch_median','mismatch_min', ...
        'n_pairs_le_1_of_1200','n_pairs_ge_90pct'});
    writetable(TrepeatRel,fullfile(outputRoot,'Mini5Pro_V4_repeat_by_rel_summary.csv'));
end
%% ================================================================
% 9. State summary (flight / hover / static)
% ================================================================
if exist('Tfile','var')
    stateNames=unique(Tfile.state,'stable');
    ST_state={}; ST_files=[]; ST_bursts=[]; ST_center_med=[];
    ST_n64=[]; ST_evm64_med=[]; ST_nq=[]; ST_evmm1_med=[]; ST_m1_min=[]; ST_allq=[];
    for ii=1:numel(stateNames)
        st=stateNames{ii};
        fsel=strcmp(Tfile.state,st);
        ST_state{end+1,1}=st; %#ok<SAGROW>
        ST_files(end+1,1)=sum(fsel); %#ok<SAGROW>
        ST_bursts(end+1,1)=sum(Tfile.n_burst(fsel)); %#ok<SAGROW>
        ST_center_med(end+1,1)=median(Tfile.center_rel_hz(fsel)); %#ok<SAGROW>
        if exist('Tsym','var')
            ssel=strcmp(Tsym.state,st);
            q64=ssel & strcmp(Tsym.modulation,'64QAM');
            qq=ssel & strcmp(Tsym.modulation,'QPSK') & (Tsym.m==1);
            ST_n64(end+1,1)=sum(q64); %#ok<SAGROW>
            if any(q64), ST_evm64_med(end+1,1)=median(Tsym.evm(q64)); else, ST_evm64_med(end+1,1)=NaN; end %#ok<SAGROW>
            ST_nq(end+1,1)=sum(qq); %#ok<SAGROW>
            if any(qq), ST_evmm1_med(end+1,1)=median(Tsym.evm(qq)); else, ST_evmm1_med(end+1,1)=NaN; end %#ok<SAGROW>
        else
            ST_n64(end+1,1)=0; ST_evm64_med(end+1,1)=NaN; %#ok<SAGROW>
            ST_nq(end+1,1)=0; ST_evmm1_med(end+1,1)=NaN; %#ok<SAGROW>
        end
        if exist('Tm1','var')
            msel=strcmp(Tm1.state,st);
            if any(msel), ST_m1_min(end+1,1)=min(Tm1.min_recurrence(msel)); else, ST_m1_min(end+1,1)=NaN; end %#ok<SAGROW>
        else
            ST_m1_min(end+1,1)=NaN; %#ok<SAGROW>
        end
        if exist('Tblock','var')
            bsel=strcmp(Tblock.state,st) & strcmp(Tblock.pattern,'ALL_QPSK');
            ST_allq(end+1,1)=sum(bsel); %#ok<SAGROW>
        else
            ST_allq(end+1,1)=0; %#ok<SAGROW>
        end
    end
    Tstate=table(ST_state,ST_files,ST_bursts,ST_center_med,ST_n64,ST_evm64_med,ST_nq,ST_evmm1_med,ST_m1_min,ST_allq, ...
        'VariableNames',{'state','files','bursts','center_rel_hz_median','n_64qam_symbols', ...
        'evm64_median','n_m1_qpsk_symbols','evm_m1_median','m1_recurrence_min','n_all_qpsk_blocks'});
    writetable(Tstate,fullfile(outputRoot,'Mini5Pro_V4_state_summary.csv'));
end
%% ================================================================
% 10. Multi-file combined constellations
% ================================================================
if PLOT_CONSTELLATION && ~isempty(all64)
    figure(5); clf; set(gcf,'Color','w','Name','Mini5Pro - All 64QAM');
    qPlot=all64; if numel(qPlot)>MAX_CONSTELLATION_POINTS, qPlot=qPlot(round(linspace(1,numel(qPlot),MAX_CONSTELLATION_POINTS))); end
    scatter(real(qPlot),imag(qPlot),1,'.');
    axis equal; grid on;
    xlim([-1.4 1.4]); ylim([-1.4 1.4]);
    xlabel('I'); ylabel('Q');
    title(sprintf('Mini5Pro V4 multi-file physical-layer 64QAM (%d points)',numel(all64)));
    saveas(gcf,fullfile(outputRoot,'Mini5Pro_V4_multifile_64QAM.png'));
    drawnow;  % Keep the constellation window open for inspection (the old version closed it).
end
if PLOT_CONSTELLATION && ~isempty(allM1)
    figure(6); clf; set(gcf,'Color','w','Name','Mini5Pro - All m1 QPSK');
    qPlot=allM1; if numel(qPlot)>MAX_CONSTELLATION_POINTS, qPlot=qPlot(round(linspace(1,numel(qPlot),MAX_CONSTELLATION_POINTS))); end
    scatter(real(qPlot),imag(qPlot),2,'.');
    axis equal; grid on;
    xlim([-1.2 1.2]); ylim([-1.2 1.2]);
    xlabel('I'); ylabel('Q');
    title(sprintf('Mini5Pro V4 multi-file m1 QPSK (%d points)',numel(allM1)));
    saveas(gcf,fullfile(outputRoot,'Mini5Pro_V4_multifile_m1_QPSK.png'));
    drawnow;  % Keep the constellation window open for inspection (the old version closed it).
end
%% ================================================================
% 11. Console summary
% ================================================================
fprintf('\n============================================================\n');
fprintf(' Mini5Pro V4 FINISHED\n');
fprintf('============================================================\n');
fprintf('Input MAT files = %d\n',numel(matFiles));
if exist('Tfile','var')
    fprintf('Files processed  = %d\n',height(Tfile));
    fprintf('Bursts total     = %d\n',sum(Tfile.n_burst));
    fprintf('Duration median  = %.3f ms\n',1e3*median(Tfile.duration_s));
end
if exist('Tsym','var')
    q64=strcmp(Tsym.modulation,'64QAM');
    qq=strcmp(Tsym.modulation,'QPSK') & (Tsym.m==1);
    fprintf('64QAM symbols    = %d\n',sum(q64));
    if any(q64), fprintf('64QAM EVM median = %.3f %%\n',100*median(Tsym.evm(q64))); end
    fprintf('m1 QPSK symbols  = %d\n',sum(qq));
    if any(qq), fprintf('m1 EVM median    = %.3f %%\n',100*median(Tsym.evm(qq))); end
end
if exist('Tm1','var')
    fprintf('m1 recurrence min= %.6f\n',min(Tm1.min_recurrence));
    fprintf('m1 exact 1.0     = %d / %d\n',sum(Tm1.min_recurrence==1),height(Tm1));
end
if exist('Tblock','var')
    fprintf('13-symbol blocks = %d\n',height(Tblock));
    fprintf('NORMAL_64 blocks = %d\n',sum(strcmp(Tblock.pattern,'NORMAL_64')));
    fprintf('ALL_QPSK blocks  = %d\n',sum(strcmp(Tblock.pattern,'ALL_QPSK')));
    fprintf('MIXED blocks     = %d\n',sum(strcmp(Tblock.pattern,'MIXED')));
end
if exist('Trepeat','var')
    fprintf('Repeat pairs     = %d\n',height(Trepeat));
    fprintf('Repeat <=0.1%%   = %d\n',sum(Trepeat.symbol_mismatch_rate<=0.001));
    fprintf('Repeat >=90%%    = %d\n',sum(Trepeat.symbol_mismatch_rate>=0.90));
end
fprintf('Output root      = %s\n',outputRoot);
fprintf('\nImportant:\n');
fprintf('  * center search is automatic; do not hard-code -20.5 MHz because some files are near -28.5 MHz.\n');
fprintf('  * 64QAM I/Q level is geometry state only, not final proven DJI 6-bit labeling.\n');
fprintf('  * same-position repeat CSV is for finding stable/filler content before descrambling/FEC is frozen.\n');
fprintf('  * old v2 per-RE phase DP is intentionally NOT used.\n');
