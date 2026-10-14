%% Air2S V1.7 - Constellation only / NO disk output (pure script)
% Based on the user's Air2S_decode.m (2026-10-12).
% Removed: CRC24A/B, LTE Gold descrambling, Turbo, payload/HEX,
%          P410/Model24/Individual30 recognition and identification scoring.
% Retained: 20 MHz band search -> 30.72 MS/s -> CP/CFO -> dual-comb ZC ->
%           CIR H1/H2 -> SFBC/Alamouti -> QPSK / 64QAM soft constellations.
% No files/folders are created. Interactive plots open only AFTER all IQ processing.
% Validation: root estimator and 2x2 SFBC inverse analytically checked;
% real Air2S IQ correctness still requires signal results (rho/EVM/plots).
clc;
clear;
close all;

%% 0. Configuration
inputPath = 'E:\DJI_figure\data\data\O3\air2s\Air2s_budon_2.4_anshi_O3\iq200ms_20260311_191519_348_000050.mat';
MAX_FILES = 0;              % 0 = all .mat files when inputPath is a folder
MANUAL_ROOTS = [];          % e.g. [85 86]; [] = auto-estimate
MANUAL_CENTER_OFFSET_HZ = NaN; % e.g. 10.5e6 relative to seg.fc, NaN = auto
PLOT_CONSTELLATION = true;  % show AFTER all PHY processing
MAX_DISPLAY_POINTS = 10000; % max points per category/burst, graphics only
% No SAVE option, no outputRoot, no mkdir, no MAT/CSV/PNG writing.

%% 1. PHY parameters (same as original Air2S_decode.m)
FS_OUT=30.72e6;
NFFT=2048;
SCS=FS_OUT/NFFT;
HALF=600;
K=(-HALF:HALF).';
KI=mod(K,NFFT)+1;
EVEN=mod(K,2)==0;
ODD=~EVEN;
ZCLEN=601;
REF_SPACING=15360;
REF_TOL=5;
REF_RHO_MIN=0.60;
ROOT_PROBE_CP_MIN=0.85;
PEAK_SEP=1500;
CIR_HALF=20;
MOD_EVM_THR=0.08;

% Original adjacent-subcarrier SFBC pairing, excluding DC.
keep=(K~=0);
km=K(keep);
fullIdx=find(keep);
aa=1:2:(numel(km)-1);
bb=aa+1;
goodPair=(km(bb)-km(aa))==1;
PAIR_I=fullIdx(aa(goodPair));
PAIR_J=fullIdx(bb(goodPair));
assert(numel(PAIR_I)==600);
SFBC_A=1/sqrt(2);

% Geometrical constellations, no payload bit labeling.
modsToTest=[4 16 64 256];
CONST=cell(numel(modsToTest),1);
for mi=1:numel(modsToTest)
    M=modsToTest(mi);
    side=round(sqrt(M));
    levels=(-(side-1):2:(side-1))/sqrt((2/3)*(M-1));
    C=zeros(M,1);
    count=0;
    for ii=1:side
        for jj=1:side
            count=count+1;
            C(count)=levels(ii)+1j*levels(jj);
        end
    end
    CONST{mi}=C;
end
idealQ=CONST{1};
ideal64=CONST{3};

%% 2. Input files
if isfolder(inputPath)
    dd=dir(fullfile(inputPath,'**','*.mat'));
    matFiles=cell(numel(dd),1);
    for ii=1:numel(dd)
        matFiles{ii}=fullfile(dd(ii).folder,dd(ii).name);
    end
else
    matFiles={inputPath};
end
if isempty(matFiles)
    error('No MAT files found: %s',inputPath);
end
if MAX_FILES>0
    matFiles=matFiles(1:min(MAX_FILES,numel(matFiles)));
end
% Cache plotted points only. No figures/graphics calls before decoding ends.
plotData={};
resultFile={};
resultStatus={};

%% 3. File loop
for fi=1:numel(matFiles)
    matFile=matFiles{fi};
    [~,baseName,~]=fileparts(matFile);
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s\n',fi,numel(matFiles),matFile);
    fprintf('============================================================\n');
    fprintf('[PHY] Reading IQ ...\n');
    try
        %% ------------------------------------------------------------
        % 4.1 Read raw IQ data, sampling frequency, and RF center.
        % -------------------------------------------------------------
        s=load(matFile,'seg');
        if ~isfield(s,'seg')
            error('Expected variable seg in MAT.');
        end
        seg=s.seg;
        Fs_raw=double(seg.fs(1));
        Fc_raw=double(seg.fc(1));
        iq0=seg.iq;
        if isstruct(iq0) && isfield(iq0,'real') && isfield(iq0,'imag')
            x_raw=double([iq0.real]).' + 1j*double([iq0.imag]).';
        elseif isreal(iq0) && ismatrix(iq0) && size(iq0,2)==2
            x_raw=double(iq0(:,1))+1j*double(iq0(:,2));
        else
            x_raw=double(iq0(:));
        end
        x_raw=x_raw(:);
        x_raw=x_raw-mean(x_raw);
        x_raw=x_raw/sqrt(mean(abs(x_raw).^2)+eps);

        fprintf('Fs_raw = %.6f MS/s, Fc_raw = %.6f GHz, samples=%d\n', ...
            Fs_raw/1e6,Fc_raw/1e9,numel(x_raw));
        fprintf('Capture duration = %.3f ms (from real sample count, not filename).\n',1e3*numel(x_raw)/Fs_raw);

        %% ------------------------------------------------------------
        % 4.2 Estimate the center of the approximately 20 MHz video channel.
        % Use native 15 kHz FFT bins (6144-point FFT at 92.16 MS/s) to select frames with clear spectral edges.
        % Combine their power spectra; do not hard-code the +10.5/+17.5 MHz offsets seen in earlier captures.
        % -------------------------------------------------------------
        if isfinite(MANUAL_CENTER_OFFSET_HZ)
            f0=MANUAL_CENTER_OFFSET_HZ;
            fprintf('Using MANUAL_CENTER_OFFSET_HZ = %+.6f MHz\n',f0/1e6);
        else
            NrawFFT=round(NFFT*Fs_raw/FS_OUT);
            nFrame=floor(numel(x_raw)/NrawFFT);
            if nFrame<8
                error('Input too short for band-center estimation.');
            end

            % Limit spectral averaging to 600 frames to avoid excessive temporary memory on long captures.
            nUse=min(nFrame,600);
            frameIDs=unique(round(linspace(1,nFrame,nUse)));
            candBin=zeros(numel(frameIDs),1);
            contrast=zeros(numel(frameIDs),1);
            Pcache=zeros(NrawFFT,numel(frameIDs),'single');

            cLo=641;
            cHi=NrawFFT-640;
            centers=(cLo:cHi).';
            for kk=1:numel(frameIDs)
                fidx=frameIDs(kk);
                xx=x_raw((fidx-1)*NrawFFT+(1:NrawFFT));
                P=abs(fftshift(fft(xx,NrawFFT))).^2;
                Pcache(:,kk)=single(P);
                cs=[0;cumsum(P(:))];
                mid=(cs(centers+601)-cs(centers-600))/1201;
                low=(cs(centers-610)-cs(centers-640))/30;
                high=(cs(centers+640)-cs(centers+610))/30;
                con=10*log10(mid./max(max(low,high),eps));
                [contrast(kk),jj]=max(con);
                candBin(kk)=centers(jj);
            end

            good=contrast>10;
            if nnz(good)<3
                [~,ord]=sort(contrast,'descend');
                good=false(size(good));
                good(ord(1:max(3,round(0.30*numel(ord)))))=true;
            end
            medc=median(candBin(good));
            good=good & abs(candBin-medc)<30;
            if nnz(good)<1
                error('Band-center estimator found no coherent frames.');
            end

            Psum=sum(double(Pcache(:,good)),2);
            cs=[0;cumsum(Psum(:))];
            mid=(cs(centers+601)-cs(centers-600))/1201;
            low=(cs(centers-610)-cs(centers-640))/30;
            high=(cs(centers+640)-cs(centers+610))/30;
            con=10*log10(mid./max(max(low,high),eps));
            [~,jj]=max(con);
            cb=centers(jj);
            f0=(cb-(NrawFFT/2+1))*Fs_raw/NrawFFT;
            fprintf('Coarse Air2S band center offset = %+.6f MHz\n',f0/1e6);
        end

        %% ------------------------------------------------------------
        % 4.3 DDC, 9.45 MHz low-pass filter, and decimation to 30.72 MS/s.
        % -------------------------------------------------------------
        fprintf('[PHY] DDC / 30.72 MS/s resampling ...\n');
        nraw=(0:numel(x_raw)-1).';
        xbb=x_raw.*exp(-1j*2*pi*f0/Fs_raw*nraw);
        [rp,rq]=rat(FS_OUT/Fs_raw,1e-12);
        if rq==round(Fs_raw/FS_OUT) && rp==1
            dec=rq;
            bLP=fir1(600,9.45e6/(Fs_raw/2),kaiser(601,9));
            yf=fftfilt(bLP,double(xbb)); % double avoids MATLAB fftfilt(single) error
            yf=[yf(301:end);zeros(300,1)];
            y=yf(1:dec:end);
        else
            y=resample(xbb,rp,rq);
        end
        y=y(:)-mean(y);
        y=y/sqrt(mean(abs(y).^2)+eps);
        clear x_raw xbb nraw;
        m=(0:numel(y)-1).';

        %% ------------------------------------------------------------
        % 4.4 Coarse CP correlation: fractional CFO and integer 15 kHz-bin correction.
        % -------------------------------------------------------------
        if numel(y)<NFFT+144+1
            error('DDC output too short.');
        end
        c=y(1:end-NFFT).*conj(y(NFFT+1:end));
        % 144-sample sliding sums: O(N) instead of repeatedly convolving.
        ccum=[0;cumsum(c)];
        cs=ccum(145:end)-ccum(1:end-144);
        ep=abs(y(1:end-NFFT)).^2;
        ecum=[0;cumsum(ep)];
        e1=ecum(145:end)-ecum(1:end-144);
        ep=abs(y(NFFT+1:end)).^2;
        ecum=[0;cumsum(ep)];
        e2=ecum(145:end)-ecum(1:end-144);
        rcp=abs(cs)./(sqrt(e1.*e2)+eps);
        use=rcp>0.90;
        if any(use)
            ff=-angle(sum(cs(use)))*FS_OUT/(2*pi*NFFT);
        else
            ff=0;
        end
        yf=y.*exp(-1j*2*pi*ff/FS_OUT*m);

        % integer 15-kHz bin offset around DC
        P=zeros(NFFT,1);
        step=NFFT*30;
        for st=1:step:(numel(yf)-NFFT+1)
            P=P+abs(fftshift(fft(yf(st:st+NFFT-1),NFFT))).^2;
        end
        bestPow=-inf;
        kc=0;
        dc=NFFT/2+1;
        for q=-40:40
            lo=dc+q-600;
            hi=dc+q+600;
            if lo>=1 && hi<=NFFT
                pp=sum(P(lo:hi));
                if pp>bestPow
                    bestPow=pp;
                    kc=q;
                end
            end
        end
        ys=yf.*exp(-1j*2*pi*(kc*SCS)/FS_OUT*m);

        %% ------------------------------------------------------------
        % 4.5 Estimate dual-comb ZC roots from second-order frequency differences at CP peaks.
        % -------------------------------------------------------------
        [~,cpLocs]=findpeaks(rcp,'MinPeakHeight',ROOT_PROBE_CP_MIN, ...
            'MinPeakDistance',1000);
        rootEven=[];
        rootOdd=[];
        rootCohEven=[];
        rootCohOdd=[];
        for kk=1:numel(cpLocs)
            pp=cpLocs(kk);
            a0=pp+144;
            if a0+NFFT-1>numel(ys), continue; end
            X=fft(ys(a0:a0+NFFT-1),NFFT);
            Xa=X(KI);

            ve=Xa(1:2:end);       % K=-600 starts even -> 601 tones
            vo=Xa(2:2:end);       % odd branch -> 600 tones

            de=ve(3:end).*ve(1:end-2).*conj(ve(2:end-1)).^2;
            do=vo(3:end).*vo(1:end-2).*conj(vo(2:end-1)).^2;
            ce=abs(sum(de))/(sum(abs(de))+eps);
            co=abs(sum(do))/(sum(abs(do))+eps);
            re=mod(-angle(sum(de))/(2*pi),1)*ZCLEN;
            ro=mod(-angle(sum(do))/(2*pi),1)*ZCLEN;
            if ce>0.90 && co>0.90
                rootEven(end+1,1)=re; %#ok<SAGROW>
                rootOdd(end+1,1)=ro; %#ok<SAGROW>
                rootCohEven(end+1,1)=ce; %#ok<SAGROW>
                rootCohOdd(end+1,1)=co; %#ok<SAGROW>
            end
        end

        if isempty(MANUAL_ROOTS)
            if isempty(rootEven)
                error('No reliable dual-comb ZC root observation found.');
            end
            ue0=round(median(rootEven));
            uo0=round(median(rootOdd));
        else
            ue0=MANUAL_ROOTS(1);
            uo0=MANUAL_ROOTS(2);
        end

        %% ------------------------------------------------------------
        % 4.6 Resolve root order and +/-1 integer-bin ambiguity using maximum normalized reference correlation.
        % -------------------------------------------------------------
        fprintf('[PHY] ZC root/order and reference correlation ...\n');
        nzc=(0:ZCLEN-1).';
        rootOrders=[ue0 uo0; uo0 ue0];
        dbList=[0 1 -1];
        bestRootRho=-inf;
        bestUE=ue0; bestUO=uo0; bestDB=0;
        bestY=ys;
        bestFp=[];
        bestRhoVec=[];
        % Unit magnitude CFO rotation leaves sample power unchanged.
        ecum=[0;cumsum(abs(ys).^2)];
        pwRef=max(ecum(NFFT+1:end)-ecum(1:end-NFFT),0);

        for di=1:numel(dbList)
            db=dbList(di);
            yb=ys.*exp(-1j*2*pi*(db*SCS)/FS_OUT*m);
            for oi=1:2
                ue=rootOrders(oi,1);
                uo=rootOrders(oi,2);
                zce=exp(-1j*pi*ue*nzc.*(nzc+1)/ZCLEN);
                zco=exp(-1j*pi*uo*nzc.*(nzc+1)/ZCLEN);
                Fp=zeros(numel(K),1);
                Fp(EVEN)=zce;
                Fp(ODD)=zco(1:nnz(ODD));

                GE=zeros(NFFT,1); GO=zeros(NFFT,1);
                GE(KI(EVEN))=Fp(EVEN);
                GO(KI(ODD))=Fp(ODD);
                tE=ifft(GE); tO=ifft(GO);
                tE=tE/(norm(tE)+eps);
                tO=tO/(norm(tO)+eps);

                % fftfilt returns conv(...,'full') truncated to input length.
                % Keep NFFT:end to match conv(...,'valid') exactly.
                cEt=fftfilt(conj(flipud(tE)),double(yb));
                cOt=fftfilt(conj(flipud(tO)),double(yb));
                cE=cEt(NFFT:end);
                cO=cOt(NFFT:end);
                rhoPair=sqrt(abs(cE).^2+abs(cO).^2)./(sqrt(pwRef)+eps);
                rr=max(rhoPair);
                if rr>bestRootRho
                    bestRootRho=rr;
                    bestUE=ue; bestUO=uo; bestDB=db;
                    bestY=yb; bestFp=Fp; bestRhoVec=rhoPair;
                end
            end
        end

        ys=bestY;
        Fp=bestFp;
        rhoRef=bestRhoVec;
        centerOffsetHz=f0+ff+(kc+bestDB)*SCS;
        rfCenterHz=Fc_raw+centerOffsetHz;

        fprintf('Detected roots = (%d,%d), center offset = %+.6f MHz\n', ...
            bestUE,bestUO,centerOffsetHz/1e6);
        fprintf('Maximum normalized reference rho = %.6f\n',max(rhoRef));

        %% ------------------------------------------------------------
        % 4.7 reference peaks / 0.5ms chains
        % -------------------------------------------------------------
        [refPks,refLocs]=findpeaks(rhoRef,'MinPeakHeight',REF_RHO_MIN, ...
            'MinPeakDistance',PEAK_SEP);
        [refLocs,ord]=sort(refLocs);
        refPks=refPks(ord);
        if isempty(refLocs)
            error('No reference peak >= %.2f.',REF_RHO_MIN);
        end

        % build chains with 15360 +/-5 samples
        chains={};
        cur=refLocs(1);
        for kk=2:numel(refLocs)
            if abs(refLocs(kk)-cur(end)-REF_SPACING)<=REF_TOL
                cur(end+1)=refLocs(kk); %#ok<SAGROW>
            else
                chains{end+1}=cur; %#ok<SAGROW>
                cur=refLocs(kk);
            end
        end
        chains{end+1}=cur;

        % Only use reference rho values belonging to valid chains for statistics.
        % V1.4 convention: a continuous chain with at least two reference symbols counts as one Air2S burst.
        % Unlike Air3S, the Air2S dual-comb reference repeats every 0.5 ms within a burst.
        % Group reference peaks into chains spaced by 15360 samples, rather than counting every peak as a burst.
        validChainIdx=find(cellfun(@numel,chains)>=2);
        N_BURSTS=numel(validChainIdx);
        if N_BURSTS<1
            error('No >=2-reference 0.5ms chain found.');
        end
        burstIDByChain=nan(numel(chains),1);
        burstStartSample=zeros(N_BURSTS,1);
        burstEndSample=zeros(N_BURSTS,1);
        burstRefCount=zeros(N_BURSTS,1);
        chainRefLoc=[];
        for bb=1:N_BURSTS
            ci=validChainIdx(bb);
            burstIDByChain(ci)=bb-1;          % Burst labels start at 00, consistent with the Air3S viewer.
            burstStartSample(bb)=chains{ci}(1);
            burstEndSample(bb)=chains{ci}(end)+NFFT-1;
            burstRefCount(bb)=numel(chains{ci});
            chainRefLoc=[chainRefLoc chains{ci}]; %#ok<AGROW>
        end

        chainRefRho=zeros(numel(chainRefLoc),1);
        for kk=1:numel(chainRefLoc)
            [~,jj]=min(abs(refLocs-chainRefLoc(kk)));
            chainRefRho(kk)=refPks(jj);
        end

        fprintf('Valid chain references = %d, rho min/median/max = %.4f / %.4f / %.4f\n', ...
            numel(chainRefRho),min(chainRefRho),median(chainRefRho),max(chainRefRho));
        fprintf('Detected Air2S bursts = %d\n',N_BURSTS);
        for bb=1:N_BURSTS
            t0ms=(burstStartSample(bb)-1)/FS_OUT*1e3;
            t1ms=(burstEndSample(bb)-1)/FS_OUT*1e3;
            fprintf('  Burst %02d : refs=%d | %.3f ~ %.3f ms\n',bb-1,burstRefCount(bb),t0ms,t1ms);
        end


        %% ------------------------------------------------------------
        % 4.8 Per-chain residual CFO correction, reference CIR and SFBC decoding.
        % -------------------------------------------------------------
        symbols=struct('chain',{},'slot',{},'m',{},'mod',{},'evmQ',{},'evm64',{}, ...
            'evmBest',{},'median_cond',{},'soft',{});
        symCount=0;

        fprintf('[PHY] SFBC decoding and QPSK/64QAM decisions ...\n');
        nt=(-CIR_HALF:CIR_HALF).';
        Aall=exp(-2j*pi*(K*nt.')/NFFT);
        Aeven=exp(-2j*pi*(K(EVEN)*nt.')/NFFT);
        Aodd=exp(-2j*pi*(K(ODD)*nt.')/NFFT);

        for ci=1:numel(chains)
            pos=chains{ci};
            if numel(pos)<2, continue; end

            yy=ys;
            % CP based residual CFO over this chain
            acc=0;
            for ri=1:numel(pos)
                bodies=pos(ri);
                if ri<numel(pos)
                    for mm=1:6
                        bodies(end+1)=pos(ri)+2208+2192*(mm-1); %#ok<SAGROW>
                    end
                end
                for bi2=1:numel(bodies)
                    b0=bodies(bi2);
                    if b0-128>=1 && b0+NFFT-1<=numel(yy)
                        acc=acc+sum(yy(b0-128:b0-1).*conj(yy(b0-128+NFFT:b0-1+NFFT)));
                    end
                end
            end
            if acc~=0
                fcp=-angle(acc)*FS_OUT/(2*pi*NFFT);
            else
                fcp=0;
            end
            lo=max(1,pos(1)-3000);
            hi=min(numel(yy),pos(end)+NFFT+3000);
            nn=(lo-1:hi-1).';
            yy(lo:hi)=yy(lo:hi).*exp(-1j*2*pi*fcp/FS_OUT*nn);

            % reference H1/H2 for all refs in chain, then remove residual phase drift
            H1raw=cell(numel(pos),1);
            H2raw=cell(numel(pos),1);
            for ri=1:numel(pos)
                p0=pos(ri);
                if p0+NFFT-1>numel(yy), continue; end
                X=fft(yy(p0:p0+NFFT-1),NFFT);
                Xa=X(KI);
                he=Xa(EVEN).*conj(Fp(EVEN));
                ho=Xa(ODD).*conj(Fp(ODD));
                ce=Aeven\he;
                co=Aodd\ho;
                H1raw{ri}=Aall*ce;
                H2raw{ri}=Aall*co;
            end

            if numel(pos)>1 && ~isempty(H1raw{1})
                ph=zeros(numel(pos),1);
                for ri=1:numel(pos)
                    ph(ri)=angle(sum(conj(H1raw{1}).*H1raw{ri}));
                end
                ph=unwrap(ph);
                pp=polyfit(double(pos(:)),ph,1);
                slope=pp(1);
                yy(lo:hi)=yy(lo:hi).*exp(-1j*slope*nn);

                % Re-estimate reference channels after phase-drift compensation.
                for ri=1:numel(pos)
                    p0=pos(ri);
                    X=fft(yy(p0:p0+NFFT-1),NFFT);
                    Xa=X(KI);
                    he=Xa(EVEN).*conj(Fp(EVEN));
                    ho=Xa(ODD).*conj(Fp(ODD));
                    H1raw{ri}=Aall*(Aeven\he);
                    H2raw{ri}=Aall*(Aodd\ho);
                end
            end

            for si=1:(numel(pos)-1)
                H1=(H1raw{si}+H1raw{si+1})/2;
                H2=(H2raw{si}+H2raw{si+1})/2;
                for mm=1:6
                    b0=pos(si)+2208+2192*(mm-1);
                    if b0+NFFT-1>numel(yy), continue; end
                    X=fft(yy(b0:b0+NFFT-1),NFFT);
                    Y=X(KI);

                    % Vectorized exact 2x2 SFBC solves on 600 adjacent pairs.
                    % AA=[a b;c d], bv=[u;v], determinant formula.
                    a=SFBC_A*H1(PAIR_I);
                    b=SFBC_A*H2(PAIR_I);
                    c2=SFBC_A*conj(H2(PAIR_J));
                    d=-SFBC_A*conj(H1(PAIR_J));
                    u=Y(PAIR_I);
                    v=conj(Y(PAIR_J));
                    detAA=a.*d-b.*c2;
                    goodDet=isfinite(detAA) & abs(detAA)>1e-14;
                    if ~all(goodDet)
                        % One unreliable 2x2 channel pair contaminates normalization/EVM.
                        % Skip this OFDM symbol instead of fabricating huge points.
                        continue;
                    end
                    s0=(d.*u-b.*v)./detAA;
                    s1=(-c2.*u+a.*v)./detAA;
                    % Spectral 2-norm condition number of each 2x2 matrix.
                    frob2=abs(a).^2+abs(b).^2+abs(c2).^2+abs(d).^2;
                    discr=sqrt(max(frob2.^2-4*abs(detAA).^2,0));
                    lmax=(frob2+discr)/2;
                    lmin=2*abs(detAA).^2./(frob2+discr+eps);
                    cnd=sqrt(lmax./max(lmin,eps));
                    st=zeros(1200,1);
                    st(1:2:end)=s0;
                    st(2:2:end)=s1;
                    zn=st/sqrt(mean(abs(st).^2)+eps);

                    metrics=zeros(numel(modsToTest),1);
                    for mdi=1:numel(modsToTest)
                        C=CONST{mdi};
                        DD=abs(zn-C.').^2;
                        [~,id1]=min(DD,[],2);
                        q=C(id1);
                        metrics(mdi)=sqrt(mean(abs(zn-q).^2)/(mean(abs(q).^2)+eps));
                    end

                    modName='OTHER';
                    for mdi=1:numel(modsToTest)
                        if metrics(mdi)<MOD_EVM_THR
                            M=modsToTest(mdi);
                            if M==4, modName='QPSK';
                            elseif M==16, modName='16QAM';
                            elseif M==64, modName='64QAM';
                            else, modName='256QAM';
                            end
                            break;
                        end
                    end

                    symCount=symCount+1;
                    symbols(symCount).chain=ci-1;
                    symbols(symCount).slot=si-1;
                    symbols(symCount).m=mm;
                    symbols(symCount).mod=modName;
                    symbols(symCount).evmQ=metrics(1);
                    symbols(symCount).evm64=metrics(3);
                    symbols(symCount).evmBest=min(metrics);
                    symbols(symCount).median_cond=median(cnd);
                    symbols(symCount).soft=single(zn);

                end
            end
        end

        if isempty(symbols)
            error('No SFBC-decoded data symbols.');
        end

        % V1.4: buffer actual recovered soft constellation points by detected burst.
        % symbols.chain is a zero-based chain index; burstIDByChain maps valid chains to Burst 00, 01, etc.
        burstQpskCell=cell(1,N_BURSTS);
        burstQamCell=cell(1,N_BURSTS);
        burstOtherCell=cell(1,N_BURSTS);
        burstQpskSymCount=zeros(1,N_BURSTS);
        burstQamSymCount=zeros(1,N_BURSTS);
        burstOtherSymCount=zeros(1,N_BURSTS);
        burst16Count=zeros(1,N_BURSTS);
        burst256Count=zeros(1,N_BURSTS);
        burstUnknownCount=zeros(1,N_BURSTS);
        burstCondCell=cell(1,N_BURSTS);
        burstQpskEvmCell=cell(1,N_BURSTS);
        burstQamEvmCell=cell(1,N_BURSTS);
        % Store both QPSK-fit and 64QAM-fit EVM for every decoded symbol; zero classified 64QAM does not imply no decoded data.
        burstAllEvmQCell=cell(1,N_BURSTS);
        burstAllEvm64Cell=cell(1,N_BURSTS);
        for si2=1:numel(symbols)
            chain1=double(symbols(si2).chain)+1;
            if chain1<1 || chain1>numel(burstIDByChain) || isnan(burstIDByChain(chain1))
                continue;
            end
            bb=burstIDByChain(chain1)+1;
            burstAllEvmQCell{bb}(end+1)=double(symbols(si2).evmQ); %#ok<AGROW>
            burstAllEvm64Cell{bb}(end+1)=double(symbols(si2).evm64); %#ok<AGROW>
            burstCondCell{bb}(end+1)=double(symbols(si2).median_cond); %#ok<AGROW>
            if strcmp(symbols(si2).mod,'QPSK')
                burstQpskCell{bb}=[burstQpskCell{bb};double(symbols(si2).soft(:))]; %#ok<AGROW>
                burstQpskSymCount(bb)=burstQpskSymCount(bb)+1;
                burstQpskEvmCell{bb}(end+1)=double(symbols(si2).evmQ); %#ok<AGROW>
            elseif strcmp(symbols(si2).mod,'64QAM')
                burstQamCell{bb}=[burstQamCell{bb};double(symbols(si2).soft(:))]; %#ok<AGROW>
                burstQamSymCount(bb)=burstQamSymCount(bb)+1;
                burstQamEvmCell{bb}(end+1)=double(symbols(si2).evm64); %#ok<AGROW>
            else
                % Not claimed as QPSK/64QAM. Keep for a separately labeled plot.
                burstOtherCell{bb}=[burstOtherCell{bb};double(symbols(si2).soft(:))]; %#ok<AGROW>
                burstOtherSymCount(bb)=burstOtherSymCount(bb)+1;
                switch symbols(si2).mod
                    case '16QAM'
                        burst16Count(bb)=burst16Count(bb)+1;
                    case '256QAM'
                        burst256Count(bb)=burst256Count(bb)+1;
                    otherwise
                        burstUnknownCount(bb)=burstUnknownCount(bb)+1;
                end
            end
        end


        % Keep display memory bounded while preserving all PHY statistics.
        for bb=1:N_BURSTS
            if numel(burstQpskCell{bb})>MAX_DISPLAY_POINTS
                iPlot=unique(round(linspace(1,numel(burstQpskCell{bb}),MAX_DISPLAY_POINTS)));
                burstQpskCell{bb}=burstQpskCell{bb}(iPlot);
            end
            if numel(burstQamCell{bb})>MAX_DISPLAY_POINTS
                iPlot=unique(round(linspace(1,numel(burstQamCell{bb}),MAX_DISPLAY_POINTS)));
                burstQamCell{bb}=burstQamCell{bb}(iPlot);
            end
            if numel(burstOtherCell{bb})>MAX_DISPLAY_POINTS
                iPlot=unique(round(linspace(1,numel(burstOtherCell{bb}),MAX_DISPLAY_POINTS)));
                burstOtherCell{bb}=burstOtherCell{bb}(iPlot);
            end
        end

        fprintf('\n[PHY] File complete: %d bursts, %d decoded OFDM data symbols.\n', N_BURSTS, numel(symbols));
        fprintf('      roots=(%d,%d), carrier=%.6f GHz, ref rho median=%.4f\n', ...
            bestUE,bestUO,rfCenterHz/1e9,median(chainRefRho));
        for bb=1:N_BURSTS
            qE=burstQpskEvmCell{bb};
            aE=burstQamEvmCell{bb};
            if isempty(qE), qPct=NaN; else, qPct=100*median(qE); end
            if isempty(aE), aPct=NaN; else, aPct=100*median(aE); end
            if isempty(burstCondCell{bb}), medCond=NaN;
            else, medCond=median(burstCondCell{bb}); end
            fprintf('      Burst %02d: QPSK=%d, 64QAM=%d, 16QAM=%d, 256QAM=%d, unknown=%d OFDM symbols\n', ...
                bb-1,burstQpskSymCount(bb),burstQamSymCount(bb),burst16Count(bb), ...
                burst256Count(bb),burstUnknownCount(bb));
            fprintf('        matched EVM (QPSK/64QAM) = %.2f%% / %.2f%%, median SFBC cond = %.2f\n', ...
                qPct,aPct,medCond);
            fprintf('        all-symbol median geometric EVM candidates = %.2f%% (QPSK), %.2f%% (64QAM)\n', ...
                100*median(burstAllEvmQCell{bb}),100*median(burstAllEvm64Cell{bb}));
        end
        p=struct();
        p.name=baseName;
        p.q=burstQpskCell;
        p.a=burstQamCell;
        p.other=burstOtherCell;
        p.qCount=burstQpskSymCount;
        p.aCount=burstQamSymCount;
        p.otherCount=burstOtherSymCount;
        p.n16=burst16Count;
        p.n256=burst256Count;
        p.nUnknown=burstUnknownCount;
        plotData{end+1}=p; %#ok<SAGROW>
        resultFile{end+1,1}=baseName; %#ok<SAGROW>
        resultStatus{end+1,1}='OK'; %#ok<SAGROW>
    catch ME
        fprintf(2,'[ERROR] %s: %s\n',baseName,ME.message);
        % Record error without stopping other files.
        resultFile{end+1,1}=baseName; %#ok<SAGROW>
        resultStatus{end+1,1}=ME.message; %#ok<SAGROW>
    end
end

fprintf('\n============================================================\n');
fprintf('Air2S PHY demodulation FINISHED (no CRC, no Turbo, no HEX).\n');
fprintf('Files processed: %d; successful: %d\n',numel(matFiles),numel(plotData));
fprintf('============================================================\n');

%% 4. Only now create the on-screen constellation figures.
% Every plot is based on actual SFBC-decoded symbol values; ideal references
% are drawn with x markers and are NOT counted as recovered points.
if PLOT_CONSTELLATION && ~isempty(plotData)
    for piFile=1:numel(plotData)
        p=plotData{piFile};
        nB=numel(p.q);
        qAll=vertcat(p.q{:});
        aAll=vertcat(p.a{:});
        otherAll=vertcat(p.other{:});
        % Upper bound the merged display too; PHY counts remain unchanged.
        if numel(qAll)>MAX_DISPLAY_POINTS
            qAll=qAll(unique(round(linspace(1,numel(qAll),MAX_DISPLAY_POINTS))));
        end
        if numel(aAll)>MAX_DISPLAY_POINTS
            aAll=aAll(unique(round(linspace(1,numel(aAll),MAX_DISPLAY_POINTS))));
        end
        if numel(otherAll)>MAX_DISPLAY_POINTS
            otherAll=otherAll(unique(round(linspace(1,numel(otherAll),MAX_DISPLAY_POINTS))));
        end
        qPages=[{qAll},p.q];
        aPages=[{aAll},p.a];
        otherPages=[{otherAll},p.other];
        labels=cell(1,nB+1);
        labels{1}=sprintf('ALL | QPSK=%d | 64QAM=%d | other=%d', ...
            sum(p.qCount),sum(p.aCount),sum(p.otherCount));
        for bb=1:nB
            labels{bb+1}=sprintf('Burst %02d | Q=%d, 64=%d, other=%d', ...
                bb-1,p.qCount(bb),p.aCount(bb),p.otherCount(bb));
        end
        % Build hidden and make visible once completed: no drawnow in PHY.
        fig=figure('Visible','off','Color','w', ...
            'Name',['Air2S Constellations - ' p.name], ...
            'NumberTitle','off','Position',[55 65 1490 625]);
        axQ=axes('Parent',fig,'Position',[0.055 0.13 0.275 0.73]);
        axA=axes('Parent',fig,'Position',[0.375 0.13 0.275 0.73]);
        axO=axes('Parent',fig,'Position',[0.695 0.13 0.275 0.73]);
        uicontrol('Parent',fig,'Style','text','Units','normalized', ...
            'Position',[0.05 0.925 0.09 0.045], ...
            'String','Select Burst:','BackgroundColor','w','FontSize',11);
        pop=uicontrol('Parent',fig,'Style','popupmenu','Units','normalized', ...
            'Position',[0.15 0.925 0.55 0.045], ...
            'String',labels,'Value',1,'FontSize',10);
        chk=uicontrol('Parent',fig,'Style','checkbox','Units','normalized', ...
            'Position',[0.735 0.925 0.25 0.045], ...
            'String','Ideal constellation (black x)','BackgroundColor','w','Value',1);
        ud=struct();
        ud.q=qPages; ud.a=aPages; ud.o=otherPages;
        ud.pop=pop; ud.chk=chk;
        ud.axQ=axQ; ud.axA=axA; ud.axO=axO;
        ud.refQ=idealQ; ud.refA=ideal64;
        set(fig,'UserData',ud);
        % Callbacks only fire when the user explicitly selects a burst.
        cb=[ ...
            'f=ancestor(gcbo,''figure''); u=get(f,''UserData''); k=get(u.pop,''Value''); ' ...
            'q=u.q{k}; a=u.a{k}; o=u.o{k}; showRef=get(u.chk,''Value'')==1; ' ...
            'cla(u.axQ); if ~isempty(q), plot(u.axQ,real(q),imag(q),''.'',''MarkerSize'',3); ' ...
            'else, text(u.axQ,0,0,''No QPSK-classified symbols'',''HorizontalAlignment'',''center''); end; ' ...
            'hold(u.axQ,''on''); if showRef, plot(u.axQ,real(u.refQ),imag(u.refQ),''kx''); end; ' ...
            'hold(u.axQ,''off''); grid(u.axQ,''on''); axis(u.axQ,''equal''); ' ...
            'xlim(u.axQ,[-1.5 1.5]); ylim(u.axQ,[-1.5 1.5]); ' ...
            'xlabel(u.axQ,''I''); ylabel(u.axQ,''Q''); title(u.axQ,''Classified QPSK''); ' ...
            'cla(u.axA); if ~isempty(a), plot(u.axA,real(a),imag(a),''.'',''MarkerSize'',3); ' ...
            'else, text(u.axA,0,0,''No 64QAM-classified symbols'',''HorizontalAlignment'',''center''); end; ' ...
            'hold(u.axA,''on''); if showRef, plot(u.axA,real(u.refA),imag(u.refA),''kx''); end; ' ...
            'hold(u.axA,''off''); grid(u.axA,''on''); axis(u.axA,''equal''); ' ...
            'xlim(u.axA,[-1.5 1.5]); ylim(u.axA,[-1.5 1.5]); ' ...
            'xlabel(u.axA,''I''); ylabel(u.axA,''Q''); title(u.axA,''Classified 64QAM''); ' ...
            'cla(u.axO); if ~isempty(o), plot(u.axO,real(o),imag(o),''.'',''MarkerSize'',3); ' ...
            'else, text(u.axO,0,0,''No other symbols'',''HorizontalAlignment'',''center''); end; ' ...
            'grid(u.axO,''on''); axis(u.axO,''equal''); ' ...
            'xlim(u.axO,[-1.5 1.5]); ylim(u.axO,[-1.5 1.5]); ' ...
            'xlabel(u.axO,''I''); ylabel(u.axO,''Q''); title(u.axO,''16/256QAM or unknown''); ' ...
        ];
        set(pop,'Callback',cb);
        set(chk,'Callback',cb);
        % Initial page is ALL bursts, drawn without using callback/gcbo.
        q=qPages{1}; a=aPages{1}; o=otherPages{1};
        if ~isempty(q)
            plot(axQ,real(q),imag(q),'.','MarkerSize',3);
        else
            text(axQ,0,0,'No QPSK-classified symbols','HorizontalAlignment','center');
        end
        hold(axQ,'on'); plot(axQ,real(idealQ),imag(idealQ),'kx'); hold(axQ,'off');
        grid(axQ,'on'); axis(axQ,'equal'); xlim(axQ,[-1.5 1.5]); ylim(axQ,[-1.5 1.5]);
        xlabel(axQ,'I'); ylabel(axQ,'Q'); title(axQ,'Classified QPSK');
        if ~isempty(a)
            plot(axA,real(a),imag(a),'.','MarkerSize',3);
        else
            text(axA,0,0,'No 64QAM-classified symbols','HorizontalAlignment','center');
        end
        hold(axA,'on'); plot(axA,real(ideal64),imag(ideal64),'kx'); hold(axA,'off');
        grid(axA,'on'); axis(axA,'equal'); xlim(axA,[-1.5 1.5]); ylim(axA,[-1.5 1.5]);
        xlabel(axA,'I'); ylabel(axA,'Q'); title(axA,'Classified 64QAM');
        if ~isempty(o)
            plot(axO,real(o),imag(o),'.','MarkerSize',3);
        else
            text(axO,0,0,'No other symbols','HorizontalAlignment','center');
        end
        grid(axO,'on'); axis(axO,'equal'); xlim(axO,[-1.5 1.5]); ylim(axO,[-1.5 1.5]);
        xlabel(axO,'I'); ylabel(axO,'Q'); title(axO,'16/256QAM or unknown');
        set(fig,'Visible','on');
    end
    fprintf('[PLOT] Opened %d constellation window(s), all remain open.\n',numel(plotData));
elseif PLOT_CONSTELLATION
    warning('No constellation can be displayed; check [ERROR] diagnostics above.');
end
