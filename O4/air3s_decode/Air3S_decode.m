%% Air3S

clc;
clear;
close all;

%% ================================================================
% 0. User configuration
% ================================================================

% Input file (V6-style configuration)
matFile = 'E:\DJI_figure\data\data\O4\air3s\Air3s_悬停\iq200ms_20260128_234612_716_000013.mat';
% Load the MAT file only once in Section 9.1 to avoid duplicated IQ arrays.
% Accept source sampling rates such as 92.16 MS/s; convert to 30.72 MS/s.

% Candidate RF channel center frequencies validated on the current dataset.
% Add other channel center frequencies here if new recordings require them.
rfCandidatesHz = [2.4205e9, 2.4405e9, 2.4605e9];

% Show constellation plots (true = enabled).
PLOT_CONSTELLATION = true;

% Display only; legacy CSV switch is retained but intentionally inactive.
SAVE_BLOCK_5700_CSV = false;

% Minimum transition agreement for high-confidence Gold pilot anchors.
GOLD_TRANSITION_MIN = 0.98;

% ZC normalized-correlation detection threshold.
ZC_RHO_MIN = 0.58;

% Complete-block F/D anchor thresholds, fixed from a 40-file evaluation.
% ser < F_SER_MAX -> F
% ser > D_SER_MIN -> D
% Values between the two thresholds are classified as U (uncertain).
F_SER_MAX = 0.20;
D_SER_MIN = 0.80;

%% ================================================================
% 1. Input file list
% ================================================================

% V6-style input: single MAT file; no result directory is created.
matFiles = {matFile};

% Display figures only; do not create output folders or write results to disk.

scriptDir = fileparts(mfilename('fullpath'));

%% ================================================================
% 2. PHY constants
% ================================================================

FS_OUT = 30.72e6;
NFFT   = 1024;
SCS    = FS_OUT/NFFT;       % 30 kHz

CP_N = 72;
CP_L = 88;

L    = 601;
ROOT = 330;

DATA_GAP = 1352;            % pFirstCP - pZC = (72+1024)+256

DC  = NFFT/2 + 1;
ACT = (DC-300):(DC+300);
K   = (-300:300).';

% With the current frequency-bin convention, odd K carries Gold QPSK pilots.
pilotIdx = find(mod(K,2)~=0);             % 300
pilotK   = K(pilotIdx);

% mixed pilot-bearing symbol: even-K, DC excluded -> 300 64QAM
mixedPayloadIdx = find(mod(K,2)==0 & K~=0);

% full 64QAM symbol: all active non-DC -> 600 64QAM
fullPayloadIdx = find(K~=0);

assert(numel(pilotIdx)==300);
assert(numel(mixedPayloadIdx)==300);
assert(numel(fullPayloadIdx)==600);

%% ================================================================
% 3. ZC reference
% ================================================================

nzc = (0:L-1).';
ZC = exp(-1j*pi*ROOT*nzc.*(nzc+1)/L);

Szc = zeros(NFFT,1);
Szc(ACT) = ZC;

tzc = ifft(ifftshift(Szc));
hZC = [tzc(end-CP_N+1:end); tzc];
EZC = sum(abs(hZC).^2);

%% ================================================================
% 4. Gold QPSK pilot: c_init = 0x34401E14
% ================================================================

GOLD_CINIT = uint32(hex2dec('34401E14'));

Nc    = 1600;
Ngold = 600;
Ntot  = Nc + Ngold + 31;

x1g = zeros(Ntot,1,'uint8');
x2g = zeros(Ntot,1,'uint8');

x1g(1)=1;

for ii=0:30
    x2g(ii+1)=uint8(bitget(GOLD_CINIT,ii+1));
end

for ii=1:(Ntot-31)
    x1g(ii+31) = bitxor(x1g(ii+3),x1g(ii));

    a1 = bitxor(x2g(ii+3),x2g(ii+2));
    a2 = bitxor(x2g(ii+1),x2g(ii));
    x2g(ii+31) = bitxor(a1,a2);
end

cGold = bitxor(x1g(Nc+1:Nc+Ngold),x2g(Nc+1:Nc+Ngold));

% The LTE-style Gray QPSK convention was selected from 24 mappings using Gold-sequence recurrence.
Rgold = ((1-2*double(cGold(1:2:end))) + ...
         1j*(1-2*double(cGold(2:2:end))))/sqrt(2);

absolute600_gold = char(double(cGold(:).')+'0');

Dgold = angle(Rgold(2:end).*conj(Rgold(1:end-1)));
DgoldClass = mod(round(Dgold/(pi/2)),4);

% Do not write Air3S_Gold_absolute600.txt; keep absolute600_gold for inspection.

%% ================================================================
% 5. Geometric 64QAM constellation (final six-bit labeling is not asserted).
% ================================================================

lev = [-7 -5 -3 -1 1 3 5 7]/sqrt(42);

C64 = zeros(64,1);
c64I = zeros(64,1);
c64Q = zeros(64,1);

ccnt=0;
for ii=1:8
    for jj=1:8
        ccnt=ccnt+1;
        C64(ccnt)=lev(ii)+1j*lev(jj);
        c64I(ccnt)=ii-1;     % Zero-based geometric level index (0..7), not a proven bit label.
        c64Q(ccnt)=jj-1;
    end
end

%% ================================================================
% 6. Optional FFF template learned from 40 files in V11.9
% ================================================================

fffTemplate=[];
fffPath=fullfile(scriptDir,'Air3S_V11_10_FFF_full5700_majority_decision.txt');

if exist(fffPath,'file')
    fffTemplate=load(fffPath);
    fffTemplate=fffTemplate(:);
    if numel(fffTemplate)~=5700
        warning('FFF template length=%d, expected5700. FFF/DDD disabled.',numel(fffTemplate));
        fffTemplate=[];
    else
        fprintf('Loaded FFF 5700-QAM template.\n');
    end
else
    fprintf('[INFO] Optional FFF template not found; FFF/DDD classification disabled (constellation decoding remains enabled).\n');
end

% Serialized layout of a 5700-QAM block:
% m2  : 1:300
% m3  : 301:900
% m4  : 901:1500
% m5  : 1501:2100
% m6  : 2101:2400
% m7  : 2401:3000
% m8  : 3001:3600
% m9  : 3601:4200
% m10 : 4201:4500
% m11 : 4501:5100
% m12 : 5101:5700
segM2  = 1:300;
segM6  = 2101:2400;
segM10 = 4201:4500;

%% ================================================================
% 7. Optional m1 transition codebook
% ================================================================

m1Codebook=[];
m1ClusterID=[];

m1CodebookPath=fullfile(scriptDir,'Air3S_V11_10_m1_transition_codebook.csv');

if exist(m1CodebookPath,'file')
    Tc=readtable(m1CodebookPath);
    m1ClusterID=Tc.cluster;
    m1Codebook=table2array(Tc(:,2:end));
    fprintf('Loaded m1 transition codebook: %d clusters.\n',size(m1Codebook,1));
end

%% ================================================================
% 8. Global result containers
% ================================================================

G_file={};
G_burst=[];
G_zc_rho=[];
G_cp_rho=[];
G_cfo_hz=[];
G_rf_center_hz=[];
G_nanchor=[];

S_file={};
S_burst=[];
S_symbol=[];
S_mod12=[];
S_ncar=[];
S_evm_raw=[];
S_evm_dd=[];

M1_file={};
M1_burst=[];
M1_block=[];
M1_rms_deg=[];
M1_cluster=[];
M1_hamming=[];
M1_abs_evm_odd=[];
M1_abs_evm_even=[];

B_file={};
B_burst=[];
B_block=[];
B_pattern={};
B_full_ser=[];
B_ser_m2=[];
B_ser_m6=[];
B_ser_m10=[];
B_nqam=[];

allPayload64=[];
allM1Qpsk=[];

%% ================================================================
% 9. Input file processing loop
% ================================================================

for fi=1:numel(matFiles)

    matFile=matFiles{fi};
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s\n',fi,numel(matFiles),matFile);
    fprintf('============================================================\n');

    fprintf('[INFO] Loading MAT...\n'); tLoad=tic;
    s=load(matFile);
    fprintf('[INFO] MAT loaded in %.2f s\n',toc(tLoad));

    % ------------------------------------------------------------
    % 9.1 Read IQ data and sample/RF frequencies
    % ------------------------------------------------------------
    if isfield(s,'seg')
        x_raw=s.seg.iq;
        Fs_raw=double(s.seg.fs(1));
        Fc_raw=double(s.seg.fc(1));
    else
        error('Current V11.11 expects seg.iq / seg.fs / seg.fc in %s',matFile);
    end

    if isreal(x_raw) && ismatrix(x_raw) && size(x_raw,2)==2
        x_raw=double(x_raw(:,1))+1j*double(x_raw(:,2));
    else
        x_raw=double(x_raw(:));
    end

    x_raw=x_raw(:);
    clear s;  % Discard the second copy of seg.IQ to reduce memory consumption.
    x_raw=x_raw-mean(x_raw);
    x_raw=x_raw./sqrt(mean(abs(x_raw).^2)+eps);
    x_raw=single(x_raw); % Retain single precision for the source IQ/STFT and RF mixing to reduce peak memory.

    [~,baseName,~]=fileparts(matFile);
    fprintf('[INFO] Samples=%d, Fs=%.3f MS/s, duration=%.1f ms\n', ...
        numel(x_raw),Fs_raw/1e6,numel(x_raw)/Fs_raw*1e3);

    % ------------------------------------------------------------
    % Render Figure 1 before the RF scan and full ZC correlation.
    % Sparse STFT: at most 1400 windows, avoiding a very large 4096-by-18000 matrix.
    % Each column is an FFT of contiguous raw IQ samples, not an aliased decimation.
    % ------------------------------------------------------------
    fprintf('[INFO] Drawing memory-bounded spectrogram...\n');
    tFig=tic;
    TF_NFFT=2048;
    TF_WIN=2048;
    TF_MAX_COLUMNS=1400;
    if numel(x_raw)<TF_WIN
        error('Not enough input samples to draw the time-frequency plot: %d',numel(x_raw));
    end
    TF_HOP=max(TF_WIN,ceil((numel(x_raw)-TF_WIN)/max(1,TF_MAX_COLUMNS-1)));
    TF_IDX=1:TF_HOP:(numel(x_raw)-TF_WIN+1);
    TF_W=single(hann(TF_WIN,'periodic'));
    P_tf=zeros(TF_NFFT,numel(TF_IDX),'single');
    for kTF=1:numel(TF_IDX)
        zz=x_raw(TF_IDX(kTF):TF_IDX(kTF)+TF_WIN-1).*TF_W;
        spec=fftshift(fft(zz,TF_NFFT));
        P_tf(:,kTF)=20*log10(abs(spec)+1e-8);
    end
    F_tf=(-TF_NFFT/2:TF_NFFT/2-1).'*(Fs_raw/TF_NFFT);
    T_tf=(TF_IDX-1+(TF_WIN-1)/2)/Fs_raw;
    if Fc_raw~=0
        fPlot_tf=(Fc_raw+F_tf)/1e6;
        yLabel_tf='RF Frequency (MHz)';
    else
        fPlot_tf=F_tf/1e6;
        yLabel_tf='Baseband Offset (MHz)';
    end
    fig1=figure(1); clf(fig1);
    set(fig1,'Color','w','Name','Figure 1 - Air3S Time-Frequency', ...
        'NumberTitle','off','Position',[80 80 1450 680]);
    axTF=axes('Parent',fig1);
    imagesc(axTF,T_tf*1e3,fPlot_tf,P_tf); axis(axTF,'xy');
    xlabel(axTF,'Time (ms)'); ylabel(axTF,yLabel_tf);
    title(axTF,sprintf('Air3S Time-Frequency | %s',baseName),'Interpreter','none');
    colormap(axTF,turbo); colorbar(axTF);
    drawnow;
    fprintf('[INFO] Figure 1 visible (%.2f s).\n',toc(tFig));
    clear P_tf spec zz;

    % ------------------------------------------------------------
    % 9.2 Candidate RF-center scan
    % ------------------------------------------------------------
    [rp,rq]=rat(FS_OUT/Fs_raw,1e-12);

    bestScore=-inf;
    bestRFHz=NaN;
    bestDDCHz=NaN;
    x=[];
    rhoZC=[];
    hMatch=conj(flipud(hZC));
    LM=numel(hMatch);

    % Evaluate several separated short windows instead of scanning entire IQ for each candidate.
    % Each window is nominally 12 ms; increase SCAN_WINDOW_MS/SCAN_COUNT if detection is weak.
    SCAN_COUNT=5;
    % A file named 200ms may contain only 20 ms; reduce overlap for short recordings.
    durationMs=1000*numel(x_raw)/Fs_raw;
    SCAN_WINDOW_MS=min(12,max(3,durationMs/SCAN_COUNT));
    scanN=min(numel(x_raw),max(round(SCAN_WINDOW_MS/1000*Fs_raw), ...
        ceil(2*LM*Fs_raw/FS_OUT)));
    scanStarts=unique(round(linspace(1,numel(x_raw)-scanN+1,SCAN_COUNT)));

    fprintf('[INFO] RF scan: %d windows of %.1f ms, %d candidates\n', ...
        numel(scanStarts),scanN/Fs_raw*1000,numel(rfCandidatesHz));
    tScan=tic;
    for ci=1:numel(rfCandidatesHz)
        fRF=rfCandidatesHz(ci);
        fDDC=fRF-Fc_raw;
        if abs(fDDC)>0.5*Fs_raw
            fprintf('  RF %.4f GHz : outside Nyquist\n',fRF/1e9);
            continue;
        end
        scanScore=-inf;
        for wi=1:numel(scanStarts)
            i1=scanStarts(wi);
            i2=i1+scanN-1;
            segScan=x_raw(i1:i2);
            nScan=(i1-1:i2-1).';
            segMix=segScan.*single(exp(-1j*2*pi*(fDDC/Fs_raw)*nScan));
            segRes=resample(segMix,rp,rq);
            % Some MATLAB releases do not support single-precision inputs to fftfilt.
            % Convert only each scan window to double; keep the large source IQ in single precision.
            segRes=double(segRes(:));
            segRes=segRes-mean(segRes);
            segRes=segRes./sqrt(mean(abs(segRes).^2)+eps);
            if numel(segRes)<LM
                continue;
            end
            cScan=fftfilt(hMatch,segRes);
            cScan=cScan(LM:end);
            eScan=movsum(abs(segRes).^2,[LM-1 0]);
            eScan=eScan(LM:end);
            rScan=abs(cScan)./(sqrt(eScan*EZC)+eps);
            scanScore=max(scanScore,max(rScan));
        end
        fprintf('  RF %.4f GHz : scan max ZC rho %.4f\n',fRF/1e9,scanScore);
        if scanScore>bestScore
            bestScore=scanScore;
            bestRFHz=fRF;
            bestDDCHz=fDDC;
        end
    end
    fprintf('[INFO] RF scan done in %.2f s.\n',toc(tScan));

    if ~isfinite(bestRFHz)
        warning('No candidate RF center inside acquisition Nyquist.');
        continue;
    end

    % Only downconvert/resample and correlate the selected RF candidate over the full capture.
    fprintf('[INFO] Full selected RF downconversion/resampling...\n');
    tFull=tic;
    nraw=(0:numel(x_raw)-1).';
    phasor=single(exp(-1j*2*pi*(bestDDCHz/Fs_raw)*nraw));
    xbb=x_raw.*phasor;
    clear phasor nraw;
    x=resample(xbb,rp,rq);
    clear xbb;
    % Use double for downstream fftfilt ZC matching to avoid single-input errors.
    % Only convert the resampled 30.72 MS/s IQ to double; the original IQ remains single.
    x=double(x(:));
    x=x-mean(x);
    x=x./sqrt(mean(abs(x).^2)+eps);
    fprintf('[INFO] Full resample %.2f s; ZC FFT matching...\n',toc(tFull));
    tMatch=tic;
    cc=fftfilt(hMatch,x);
    cc=cc(LM:end); % Trim the same-length convolution to match conv(...,'valid').
    ee=movsum(abs(x).^2,[LM-1 0]);
    ee=ee(LM:end);
    rhoZC=abs(cc)./(sqrt(ee*EZC)+eps);
    clear cc ee;
    bestScoreFull=max(rhoZC);
    fprintf('[INFO] Full ZC max rho=%.4f (%.2f s)\n',bestScoreFull,toc(tMatch));

    if bestScoreFull<ZC_RHO_MIN
        warning('No reliable Air3S ZC in %s; full rho %.4f (scan %.4f)', ...
            baseName,bestScoreFull,bestScore);
        title(axTF,sprintf('Air3S: No ZC found (max rho %.3f) | %s',bestScoreFull,baseName), ...
            'Interpreter','none');
        drawnow;
        continue;
    end

    fprintf('  Selected RF %.6f GHz, DDC %+.3f MHz\n',bestRFHz/1e9,bestDDCHz/1e6);

    [pksZC,locsZC]=findpeaks(rhoZC, ...
        'MinPeakHeight',ZC_RHO_MIN, ...
        'MinPeakDistance',round(1.5e-3*FS_OUT));

    fprintf('  Detected bursts = %d\n',numel(locsZC));

    % Figure 1 is already drawn; overlay ZC positions without recalculating the STFT.
    hold(axTF,'on');
    yl_tf=ylim(axTF);
    for kk=1:numel(locsZC)
        tBurst_ms=(locsZC(kk)-1)/FS_OUT*1e3;
        plot(axTF,[tBurst_ms tBurst_ms],yl_tf,'w--','LineWidth',0.7);
        % Skip per-burst labels when there are many bursts to avoid slow, crowded rendering.
        if numel(locsZC)<=30
            text(axTF,tBurst_ms,yl_tf(2)-0.04*diff(yl_tf),sprintf(' B%02d',kk-1), ...
                'Color','w','FontWeight','bold','VerticalAlignment','top');
        end
    end
    hold(axTF,'off');
    drawnow;

    % Figure 2: per-burst constellation buffers
    N_BURSTS=numel(locsZC);
    burstQpskCell=cell(1,N_BURSTS);
    burstQamCell=cell(1,N_BURSTS);
    burstEvmUI=nan(1,N_BURSTS);
    burstStatus=cell(1,N_BURSTS);
    for kk=1:N_BURSTS
        burstStatus{kk}='ZC detected; waiting for CP and Gold verification';
    end
    if PLOT_CONSTELLATION && N_BURSTS>0
        fig2=figure(2); clf(fig2);
        set(fig2,'Color','w','Name','Figure 2 - Air3S Constellation Browser', ...
            'NumberTitle','off','Position',[70 70 1250 620]);
        text(0.1,0.6,sprintf('Detected %d bursts. Recovering QPSK / 64QAM...',N_BURSTS), ...
            'FontSize',16,'Interpreter','none');
        axis off; drawnow;
    end

    % 9.3 Evaluate CP72 timing only near detected bursts to avoid full-length temporary arrays.

    % file-local payload
    filePayload64=[];
    fileM1Qpsk=[];

    fileBlockDecisions={};
    fileBlockRows=[];

    % ============================================================
    % 10. Burst loop
    % ============================================================

    for bi=1:numel(locsZC)

        zcMF=locsZC(bi);

        % --------------------------------------------------------
        % 10.1 Refine ZC timing using cyclic-prefix correlation
        % --------------------------------------------------------
        lo=max(1,zcMF-500);
        hi=min(numel(x)-NFFT-CP_N+1,zcMF+100);
        if hi<lo
            burstStatus{bi}='CP search out of signal range';
            continue;
        end
        cpA=x(lo:hi+CP_N-1);
        cpB=x(lo+NFFT:hi+NFFT+CP_N-1);
        cpNum=conv(cpA.*conj(cpB),ones(CP_N,1),'valid');
        cpEa=conv(abs(cpA).^2,ones(CP_N,1),'valid');
        cpEb=conv(abs(cpB).^2,ones(CP_N,1),'valid');
        cpRho=abs(cpNum)./(sqrt(cpEa.*cpEb)+eps);
        [rhoCp,jj]=max(cpRho);
        pZC=lo+jj-1;

        if rhoCp<0.70
            burstStatus{bi}=sprintf('CP rejected: rho %.3f < 0.70',rhoCp);
            continue;
        end

        pFirst=pZC+DATA_GAP;

        % Allow up to 48 data symbols for every burst.
        % The old code limited the last burst to seven and could discard complete blocks.
        % Actual symbol availability is still checked against capture boundaries.
        maxData=48;

        pList=zeros(maxData,1);
        cpList=zeros(maxData,1);

        pList(1)=pFirst;

        for si=1:maxData

            cpList(si)=CP_N;
            if mod(si-1,14)==0
                cpList(si)=CP_L;
            end

            if si>1
                pList(si)=pList(si-1)+NFFT+cpList(si-1);
            end
        end

        % --------------------------------------------------------
        % 10.2 all-CP CFO
        % --------------------------------------------------------
        if pZC+NFFT+CP_N-1>numel(x)
            continue;
        end

        cSum=sum(conj(x(pZC:pZC+CP_N-1)).* ...
                 x(pZC+NFFT:pZC+NFFT+CP_N-1));

        for si=1:maxData
            pp=pList(si);
            CC=cpList(si);

            if pp+NFFT+CC-1<=numel(x)
                cSum=cSum+sum(conj(x(pp:pp+CC-1)).* ...
                              x(pp+NFFT:pp+NFFT+CC-1));
            end
        end

        cfoAll=angle(cSum)*FS_OUT/(2*pi*NFFT);

        % Compensate only the current burst region; never rotate the entire 200ms recording for each burst.
        xStart=max(1,pZC-16);
        xStop=min(numel(x),pList(end)+NFFT+max(cpList)+16);
        nnLocal=(xStart:xStop).';
        xf=x(xStart:xStop).*exp(-1j*2*pi*cfoAll/FS_OUT*(nnLocal-(pZC-1)));

        % --------------------------------------------------------
        % 10.3 ZC channel estimate for Gold-anchor selection
        % --------------------------------------------------------
        aZ=pZC+CP_N;

        if aZ+NFFT-1>xStop
            burstStatus{bi}='Insufficient ZC samples';
            continue;
        end

        aZLocal=aZ-xStart+1;
        YzFull=fftshift(fft(xf(aZLocal:aZLocal+NFFT-1),NFFT));
        Yz=YzFull(ACT);

        Hls=Yz./(ZC+eps);

        % Smooth the channel estimate using a short channel impulse response (CIR).
        Hfull=zeros(NFFT,1);
        Hfull(ACT)=Hls;

        ht=ifft(ifftshift(Hfull));

        hkeep=zeros(NFFT,1);
        hkeep(1:64)=ht(1:64);
        hkeep(end-7:end)=ht(end-7:end);

        HzcFull=fftshift(fft(hkeep,NFFT));
        Hzc=HzcFull(ACT);

        % --------------------------------------------------------
        % 10.4 FFT of all available data symbols
        % --------------------------------------------------------
        Ycell=cell(maxData,1);
        tSym=nan(maxData,1);

        for si=1:maxData

            aa=pList(si)+cpList(si);

            if aa+NFFT-1<=xStop

                aaLocal=aa-xStart+1;
                yy=fftshift(fft(xf(aaLocal:aaLocal+NFFT-1),NFFT));
                Ycell{si}=yy(ACT);
                tSym(si)=aa+NFFT/2;
            end
        end

        % --------------------------------------------------------
        % 10.5 m2/m6/m10 Gold anchor -> Hdata
        % --------------------------------------------------------
        anchorSi=[];
        anchorT=[];
        anchorH={};
        anchorTransition=[];
        bestGoldPct=0;
        nGoldTried=0;

        for si=1:maxData

            mm12=mod(si-1,12)+1;

            if ~any(mm12==[2 6 10])
                continue;
            end

            if isempty(Ycell{si})
                continue;
            end

            Y=Ycell{si};

            % Check pilot phase transitions with a coarse ZC-based equalizer, without known pilot labels.
            qZc=Y(pilotIdx).*conj(Hzc(pilotIdx))./ ...
                (abs(Hzc(pilotIdx)).^2 + 0.02*mean(abs(Hzc).^2));

            dphi=angle(qZc(2:end).*conj(qZc(1:end-1)));
            dClass=mod(round(dphi/(pi/2)),4);

            transitionPct=mean(dClass==DgoldClass);
            nGoldTried=nGoldTried+1;
            bestGoldPct=max(bestGoldPct,transitionPct);

            if transitionPct<GOLD_TRANSITION_MIN
                continue;
            end

            % Coherent channel estimation from known Gold QPSK pilots.
            Hp=Y(pilotIdx)./(Rgold+eps);

            hpPhase=unwrap(angle(Hp));
            hpLogMag=log(abs(Hp)+eps);

            Hall=exp(interp1(pilotK,hpLogMag,K,'linear','extrap') + ...
                     1j*interp1(pilotK,hpPhase,K,'linear','extrap'));

            anchorSi(end+1,1)=si; %#ok<SAGROW>
            anchorT(end+1,1)=tSym(si); %#ok<SAGROW>
            anchorH{end+1,1}=Hall; %#ok<SAGROW>
            anchorTransition(end+1,1)=transitionPct; %#ok<SAGROW>
        end

        % If too few Gold anchors remain, still show a coarse-equalized diagnostic constellation.
        % Diagnostic constellations are NOT counted as validated decoding or 5700-QAM blocks.
        if numel(anchorSi)<2
            if numel(anchorSi)==1
                % One Gold anchor may support a preview but cannot validate time-domain interpolation.
                Hrough=anchorH{1};
                roughSource='single Gold anchor';
            else
                Hrough=Hzc;
                roughSource='ZC channel';
            end
            burstStatus{bi}=sprintf('Diagnostic: %s | Gold anchors=%d | best match %.1f%% (%d tried)', ...
                roughSource,numel(anchorSi),100*bestGoldPct,nGoldTried);
            if PLOT_CONSTELLATION
                for si=1:maxData
                    if isempty(Ycell{si}), continue; end
                    mm12=mod(si-1,12)+1;
                    Xrough=Ycell{si}./(Hrough+eps);
                    if mm12==1
                        zq=Xrough(K~=0);
                        zq=zq./sqrt(mean(abs(zq).^2)+eps);
                        burstQpskCell{bi}=[burstQpskCell{bi};zq]; %#ok<AGROW>
                    else
                        if any(mm12==[2 6 10])
                            ix=mixedPayloadIdx;
                        else
                            ix=fullPayloadIdx;
                        end
                        za=Xrough(ix);
                        za=za./sqrt(mean(abs(za).^2)+eps);
                        burstQamCell{bi}=[burstQamCell{bi};za]; %#ok<AGROW>
                    end
                end
            end
            if mod(bi,10)==0 || bi==1
                fprintf('[BURST %d/%d] %s\n',bi,N_BURSTS,burstStatus{bi});
                drawnow limitrate;
            end
            continue;
        end
        burstStatus{bi}=sprintf('Gold verified | anchors=%d | best match %.1f%%', ...
            numel(anchorSi),100*bestGoldPct);

        % burst summary
        G_file{end+1,1}=baseName; %#ok<SAGROW>
        G_burst(end+1,1)=bi-1; %#ok<SAGROW>
        G_zc_rho(end+1,1)=pksZC(bi); %#ok<SAGROW>
        G_cp_rho(end+1,1)=rhoCp; %#ok<SAGROW>
        G_cfo_hz(end+1,1)=cfoAll; %#ok<SAGROW>
        G_rf_center_hz(end+1,1)=bestRFHz; %#ok<SAGROW>
        G_nanchor(end+1,1)=numel(anchorSi); %#ok<SAGROW>

        % --------------------------------------------------------
        % 10.6 m1 absolute QPSK + transition cluster
        % Use the immediately following m2 Gold anchor, not Q4.
        % --------------------------------------------------------
        for bl=0:3

            s1=bl*12+1;
            s2=s1+1;

            if s2>maxData || isempty(Ycell{s1})
                continue;
            end

            ai=find(anchorSi==s2,1);

            if isempty(ai)
                continue;
            end

            Hm2=anchorH{ai};
            Xm1=Ycell{s1}./(Hm2+eps);

            % Check the odd- and even-subcarrier QPSK components separately.
            xo=Xm1(pilotIdx);
            xe=Xm1(mixedPayloadIdx);

            xo=xo./sqrt(mean(abs(xo).^2)+eps);
            xe=xe./sqrt(mean(abs(xe).^2)+eps);

            % Use nearest ideal QPSK points only for EVM computation.
            Qo=(sign(real(xo))+1j*sign(imag(xo)))/sqrt(2);
            Qe=(sign(real(xe))+1j*sign(imag(xe)))/sqrt(2);

            evOdd=sqrt(mean(abs(xo-Qo).^2)/(mean(abs(Qo).^2)+eps));
            evEven=sqrt(mean(abs(xe-Qe).^2)/(mean(abs(Qe).^2)+eps));

            fileM1Qpsk=[fileM1Qpsk;xo;xe]; %#ok<AGROW>
            burstQpskCell{bi}=[burstQpskCell{bi};xo;xe]; %#ok<AGROW>
            allM1Qpsk=[allM1Qpsk;xo;xe]; %#ok<AGROW>

            % Codebook still uses 299 transitions across the odd-K positions.
            d=angle(xo(2:end).*conj(xo(1:end-1)));
            cls=mod(round(d/(pi/2)),4);
            res=angle(exp(1j*(d-cls*pi/2)));
            rmsDeg=sqrt(mean(res.^2))*180/pi;

            bestCluster=NaN;
            bestHam=NaN;

            if ~isempty(m1Codebook)

                hd=mean(m1Codebook~=double(cls(:).'),2);
                [bestHam,jj2]=min(hd);
                bestCluster=m1ClusterID(jj2);
            end

            M1_file{end+1,1}=baseName; %#ok<SAGROW>
            M1_burst(end+1,1)=bi-1; %#ok<SAGROW>
            M1_block(end+1,1)=bl; %#ok<SAGROW>
            M1_rms_deg(end+1,1)=rmsDeg; %#ok<SAGROW>
            M1_cluster(end+1,1)=bestCluster; %#ok<SAGROW>
            M1_hamming(end+1,1)=bestHam; %#ok<SAGROW>
            M1_abs_evm_odd(end+1,1)=evOdd; %#ok<SAGROW>
            M1_abs_evm_even(end+1,1)=evEven; %#ok<SAGROW>
        end

        % --------------------------------------------------------
        % 10.7 Interpolate Hdata in time and frequency to recover 64QAM in m2 through m12.
        % --------------------------------------------------------
        decCell=cell(maxData,1);

        for si=1:maxData

            mm12=mod(si-1,12)+1;

            if mm12==1 || isempty(Ycell{si})
                continue;
            end

            tt=tSym(si);

            % Interpolate between channel anchors in the time direction.
            if tt<=anchorT(1)

                Hall=anchorH{1};

            elseif tt>=anchorT(end)

                Hall=anchorH{end};

            else

                j1=find(anchorT>=tt,1,'first');
                j0=j1-1;

                ww=(tt-anchorT(j0))/(anchorT(j1)-anchorT(j0));

                % Anchor phases are sufficiently continuous to allow direct complex interpolation.
                Hall=(1-ww)*anchorH{j0}+ww*anchorH{j1};
            end

            if any(mm12==[2 6 10])
                idx=mixedPayloadIdx;     % 300 even-K 64QAM
            else
                idx=fullPayloadIdx;      % 600 non-DC 64QAM
            end

            X=Ycell{si}(idx)./(Hall(idx)+eps);

            % Normalize each OFDM data symbol independently.
            z=X./sqrt(mean(abs(X).^2)+eps);

            % Initial nearest-neighbor geometric 64QAM decisions.
            DD=abs(z-C64.').^2;
            [~,id]=min(DD,[],2);
            s64=C64(id);

            evRaw=sqrt(mean(abs(z-s64).^2)/(mean(abs(s64).^2)+eps));

            % Two rounds of decision-directed complex-scalar refinement.
            for it=1:2

                alpha=sum(conj(s64).*z)/(sum(abs(s64).^2)+eps);
                z=z/(alpha+eps);

                DD=abs(z-C64.').^2;
                [~,id]=min(DD,[],2);
                s64=C64(id);
            end

            evDD=sqrt(mean(abs(z-s64).^2)/(mean(abs(s64).^2)+eps));

            decCell{si}=id(:)-1;       % 0-based 0..63

            filePayload64=[filePayload64;z]; %#ok<AGROW>
            burstQamCell{bi}=[burstQamCell{bi};z]; %#ok<AGROW>
            allPayload64=[allPayload64;z]; %#ok<AGROW>

            S_file{end+1,1}=baseName; %#ok<SAGROW>
            S_burst(end+1,1)=bi-1; %#ok<SAGROW>
            S_symbol(end+1,1)=si; %#ok<SAGROW>
            S_mod12(end+1,1)=mm12; %#ok<SAGROW>
            S_ncar(end+1,1)=numel(idx); %#ok<SAGROW>
            S_evm_raw(end+1,1)=evRaw; %#ok<SAGROW>
            S_evm_dd(end+1,1)=evDD; %#ok<SAGROW>
        end

        % --------------------------------------------------------
        % 10.8 Complete 12-symbol block: m2 to m12 contain 5700 QAM symbols.
        % --------------------------------------------------------
        for bl=0:3

            ss=(bl*12+2):(bl*12+12);

            if max(ss)>maxData
                continue;
            end

            good=true;

            for kk=1:numel(ss)
                if isempty(decCell{ss(kk)})
                    good=false;
                    break;
                end
            end

            if ~good
                continue;
            end

            vv=[];

            for kk=1:numel(ss)
                vv=[vv;decCell{ss(kk)}]; %#ok<AGROW>
            end

            if numel(vv)~=5700
                continue;
            end

            pattern='UNCLASSIFIED';
            fullSer=NaN;
            ser2=NaN;
            ser6=NaN;
            ser10=NaN;

            if ~isempty(fffTemplate)

                fullSer=mean(double(vv)~=double(fffTemplate));

                ser2=mean(double(vv(segM2)) ~=double(fffTemplate(segM2)));
                ser6=mean(double(vv(segM6)) ~=double(fffTemplate(segM6)));
                ser10=mean(double(vv(segM10))~=double(fffTemplate(segM10)));

                ssAnchor=[ser2 ser6 ser10];
                pp='UUU';

                for jj3=1:3
                    if ssAnchor(jj3)<F_SER_MAX
                        pp(jj3)='F';
                    elseif ssAnchor(jj3)>D_SER_MIN
                        pp(jj3)='D';
                    else
                        pp(jj3)='U';
                    end
                end

                pattern=pp;
            end

            B_file{end+1,1}=baseName; %#ok<SAGROW>
            B_burst(end+1,1)=bi-1; %#ok<SAGROW>
            B_block(end+1,1)=bl; %#ok<SAGROW>
            B_pattern{end+1,1}=pattern; %#ok<SAGROW>
            B_full_ser(end+1,1)=fullSer; %#ok<SAGROW>
            B_ser_m2(end+1,1)=ser2; %#ok<SAGROW>
            B_ser_m6(end+1,1)=ser6; %#ok<SAGROW>
            B_ser_m10(end+1,1)=ser10; %#ok<SAGROW>
            B_nqam(end+1,1)=5700; %#ok<SAGROW>

            fileBlockDecisions{end+1,1}=vv; %#ok<SAGROW>
            fileBlockRows(end+1,1)=numel(B_file); %#ok<SAGROW>
        end
        if mod(bi,10)==0 || bi==1 || bi==N_BURSTS
            fprintf('[BURST %d/%d] %s, QPSK=%d, 64QAM=%d\n', ...
                bi,N_BURSTS,burstStatus{bi}, ...
                numel(burstQpskCell{bi}),numel(burstQamCell{bi}));
            drawnow limitrate;
        end
    end

    %% ============================================================
    % 11. Current-file results
    % ============================================================

    % Do not export 5700-QAM decision-ID CSV files.

    % Do not export per-file block-summary CSV files.

    % ============================================================
    % Figure 2 initially selects the first Gold-verified burst that contains 64QAM.
    % If no Gold-verified burst exists, show the coarse-equalized diagnostic constellation marked NOT VERIFIED.
    % ============================================================
    if PLOT_CONSTELLATION && N_BURSTS>0
        figure(2); clf; fig2=gcf;
        set(fig2,'Color','w','Name','Figure 2 - Air3S Burst Constellation Browser', ...
            'NumberTitle','off','Position',[70 70 1250 620]);
        axQ=axes('Parent',fig2,'Position',[0.07 0.14 0.40 0.72]);
        axA=axes('Parent',fig2,'Position',[0.56 0.14 0.40 0.72]);
        popupLabels=cell(1,N_BURSTS);
        for kk=1:N_BURSTS
            popupLabels{kk}=sprintf('Burst %02d | Q=%d 64QAM=%d', ...
                kk-1,numel(burstQpskCell{kk}),numel(burstQamCell{kk}));
        end
        uicontrol('Parent',fig2,'Style','text','Units','normalized', ...
            'Position',[0.07 0.92 0.12 0.05],'String','Select Burst:', ...
            'BackgroundColor','w','FontSize',12,'FontWeight','bold', ...
            'HorizontalAlignment','left');
        popBurst=uicontrol('Parent',fig2,'Style','popupmenu','Units','normalized', ...
            'Position',[0.19 0.925 0.34 0.05],'String',popupLabels,'Value',1);
        chkIdeal=uicontrol('Parent',fig2,'Style','checkbox','Units','normalized', ...
            'Position',[0.56 0.925 0.30 0.05],'String','Overlay ideal 64QAM (black x)', ...
            'BackgroundColor','w','FontSize',11,'Value',1);
        uiStatus=uicontrol('Parent',fig2,'Style','text','Units','normalized', ...
            'Position',[0.07 0.03 0.89 0.06],'String','', ...
            'BackgroundColor','w','FontSize',11,'HorizontalAlignment','left');
        ud=struct('axQ',axQ,'axA',axA,'Q',{burstQpskCell}, ...
            'A',{burstQamCell},'status',{burstStatus},'qamIdeal',C64, ...
            'chkIdeal',chkIdeal,'popBurst',popBurst,'uiStatus',uiStatus);
        set(fig2,'UserData',ud);
        set(popBurst,'Callback',@(~,~) air3sRefreshConstellation(fig2));
        set(chkIdeal,'Callback',@(~,~) air3sRefreshConstellation(fig2));
        goodGold=cellfun(@(t) ~isempty(strfind(t,'Gold verified')),burstStatus);
        firstGood=find(~cellfun(@isempty,burstQamCell) & goodGold,1,'first');
        if isempty(firstGood)
            firstGood=find(~cellfun(@isempty,burstQamCell),1,'first');
        end
        if isempty(firstGood), firstGood=1; end
        set(popBurst,'Value',firstGood);
        air3sRefreshConstellation(fig2);
        drawnow;
        fprintf('[INFO] Figure 2 displayed; selected Burst %02d.\n',firstGood-1);
    end
end

%% ================================================================
% 12. Global summary tables (in memory only)
% ================================================================

if ~isempty(G_burst)

    Tburst=table(G_file,G_burst,G_rf_center_hz,G_zc_rho,G_cp_rho,G_cfo_hz,G_nanchor, ...
        'VariableNames',{'file','burst','rf_center_hz','zc_rho','cp_rho', ...
                         'cfo_hz','n_gold_anchor'});

    % Do not save CSV files.
end

if ~isempty(S_symbol)

    Tsym=table(S_file,S_burst,S_symbol,S_mod12,S_ncar,S_evm_raw,S_evm_dd, ...
        'VariableNames',{'file','burst','symbol','mod12','ncar','evm_raw','evm_dd'});

    % Do not save CSV files.
end

if ~isempty(M1_block)

    Tm1=table(M1_file,M1_burst,M1_block,M1_rms_deg,M1_cluster,M1_hamming, ...
              M1_abs_evm_odd,M1_abs_evm_even, ...
        'VariableNames',{'file','burst','block','transition_rms_deg','cluster', ...
                         'cluster_hamming','evm_odd','evm_even'});

    % Do not save CSV files.
end

if ~isempty(B_block)

    Tblock=table(B_file,B_burst,B_block,B_pattern,B_full_ser, ...
                 B_ser_m2,B_ser_m6,B_ser_m10,B_nqam, ...
        'VariableNames',{'file','burst','block','pattern','full_ser_to_FFF', ...
                         'ser_m2','ser_m6','ser_m10','nqam'});

    % Do not save CSV files.
end

% The combined Figure 3/4 was removed; keep the RF spectrogram and per-burst constellation viewer.

%% ================================================================
% 14. Console summary
% ================================================================

fprintf('\n============================================================\n');
fprintf(' Air3S V11.12 FIXED FINISHED\n');
fprintf('============================================================\n');
fprintf('Input MAT files = %d\n',numel(matFiles));

if exist('Tburst','var')
    fprintf('Valid bursts    = %d\n',height(Tburst));
    fprintf('ZC rho median   = %.4f\n',median(Tburst.zc_rho));
    fprintf('CP rho median   = %.6f\n',median(Tburst.cp_rho));
end

if exist('Tsym','var')
    fprintf('Payload symbols = %d\n',height(Tsym));
    fprintf('64QAM EVM median= %.3f %%\n',100*median(Tsym.evm_dd));
end

if exist('Tblock','var')

    fprintf('Complete blocks = %d\n',height(Tblock));

    fprintf('FFF = %d\n',sum(strcmp(Tblock.pattern,'FFF')));
    fprintf('DDD = %d\n',sum(strcmp(Tblock.pattern,'DDD')));

    mixN=sum(~strcmp(Tblock.pattern,'FFF') & ~strcmp(Tblock.pattern,'DDD'));
    fprintf('MIX/U = %d\n',mixN);
end

if exist('Tm1','var')
    fprintf('m1 observations = %d\n',height(Tm1));
    fprintf('m1 odd EVM med  = %.3f %%\n',100*median(Tm1.evm_odd));
    fprintf('m1 even EVM med = %.3f %%\n',100*median(Tm1.evm_even));
end

fprintf('Gold c_init     = 0x34401E14\n');
% No output directory is created.

fprintf('\nImportant:\n');
fprintf('  * 64QAM decisionID is geometry state 0..63, not final proven coded-bit labeling.\n');
fprintf('  * FFF/DDD is engineering classification; FFF is currently an idle/filler candidate.\n');
fprintf('  * Do not restore old 342-sample ZC shift.\n');
fprintf('  * Do not use old FIXED600 or QPSK x^4 quadrant recovery.\n');


%% ================================================================
% Local helper function (MATLAB R2016b or newer supports local functions in scripts).
% ================================================================
function air3sRefreshConstellation(fig2)
    if ~isgraphics(fig2), return; end
    ud=get(fig2,'UserData');
    k=get(ud.popBurst,'Value');
    k=max(1,min(numel(ud.Q),round(k)));
    q=ud.Q{k};
    a=ud.A{k};
    status=ud.status{k};
    isVerified=~isempty(strfind(status,'Gold verified')); %#ok<STREMP>

    cla(ud.axQ);
    if ~isempty(q)
        plot(ud.axQ,real(q),imag(q),'.','MarkerSize',4);
    end
    grid(ud.axQ,'on'); axis(ud.axQ,'equal');
    xlim(ud.axQ,[-1.5 1.5]); ylim(ud.axQ,[-1.5 1.5]);
    xlabel(ud.axQ,'I'); ylabel(ud.axQ,'Q');
    if isVerified
        title(ud.axQ,sprintf('Burst %02d: m1 QPSK (Gold verified)',k-1));
    else
        title(ud.axQ,sprintf('Burst %02d: coarse QPSK (NOT VERIFIED)',k-1));
    end

    cla(ud.axA);
    if ~isempty(a)
        plot(ud.axA,real(a),imag(a),'.','MarkerSize',3);
        hold(ud.axA,'on');
    end
    if get(ud.chkIdeal,'Value')==1
        plot(ud.axA,real(ud.qamIdeal),imag(ud.qamIdeal),'kx', ...
            'MarkerSize',6,'LineWidth',1);
    end
    hold(ud.axA,'off');
    grid(ud.axA,'on'); axis(ud.axA,'equal');
    xlim(ud.axA,[-1.5 1.5]); ylim(ud.axA,[-1.5 1.5]);
    xlabel(ud.axA,'I'); ylabel(ud.axA,'Q');
    if isVerified
        title(ud.axA,sprintf('Burst %02d: recovered 64QAM',k-1));
    else
        title(ud.axA,sprintf('Burst %02d: coarse 64QAM (NOT VERIFIED)',k-1));
    end
    set(ud.uiStatus,'String',sprintf('Burst %02d: %s | m1=%d points, 64QAM=%d points', ...
        k-1,status,numel(q),numel(a)));
end
