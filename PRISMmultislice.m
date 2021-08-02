function [emdSTEM,psiRealplot,psiFFTplot] = PRISMmultislice(emdSTEM,fbase)
tic
% Mutlislice STEM simulation to compare with PRISM.

flagPlotProbe = true;

indFP = 1;%2-1 + 3;  % Do a single frozen phonon each time
emdSTEM.MULTIprobeDefocusArray = 0;
emdSTEM.MULTIprobeSemiangleArray = 0.48/1000;%27.5/1000;  % Rads
compFreq = 0.01;  % How often to output status to console in completion fraction units

% flagPlot = 1*1;
flagPlot = true;
intRange = [0 1e-2];
intRangeFFT = [0 0.01];
flagKeep4Doutput = true;

% Probe positions
dxy = 0.25;%emdSTEM.cellDim(1) / 120;
% xR = [0.5 0.5]*emdSTEM.cellDim(1);
% yR = [0.5 0.5]*emdSTEM.cellDim(2);
% xR = [0.6 0.64]*emdSTEM.cellDim(1);
% yR = [0.6 0.64]*emdSTEM.cellDim(2);
xR = [0 1]*emdSTEM.cellDim(1);
yR = [0 1]*emdSTEM.cellDim(2);
emdSTEM.MULTIxp = (xR(1)+dxy/2):dxy:(xR(2)-dxy/2);
emdSTEM.MULTIyp = (yR(1)+dxy/2):dxy:(yR(2)-dxy/2);
% [length(emdSTEM.MULTIxp) length(emdSTEM.MULTIyp)]
emdSTEM.MULTIxp = emdSTEM.cellDim(1) * 0.5;
emdSTEM.MULTIyp = emdSTEM.cellDim(2) * 0.5;


emdSTEM.E0 = 300e3;  % Microscope voltage in volts
% Calculate wavelength and electron interaction parameter
m = 9.109383*10^-31;
e = 1.602177*10^-19;
c =  299792458;
h = 6.62607*10^-34;
emdSTEM.lambda = h/sqrt(2*m*e*emdSTEM.E0) ...
    /sqrt(1 + e*emdSTEM.E0/2/m/c^2) * 10^10; % wavelength in A
emdSTEM.sigma = (2*pi/emdSTEM.lambda/emdSTEM.E0) ...
    *(m*c^2+e*emdSTEM.E0)/(2*m*c^2+e*emdSTEM.E0);

% Detector coords
dr = 1 / 1000;
qx = makeFourierCoords(emdSTEM.imageSize(1),emdSTEM.pixelSize(1));
qy = makeFourierCoords(emdSTEM.imageSize(2),emdSTEM.pixelSize(2));
emdSTEM.qMax = min(max(abs(qx)),max(abs(qy)))/2;
alphaMax = emdSTEM.qMax * emdSTEM.lambda;
emdSTEM.MULTIdetectorAngles = (dr/2):dr:(alphaMax-dr/2);
[emdSTEM.qya,emdSTEM.qxa] = meshgrid(qy,qx);
q2 = emdSTEM.qxa.^2 + emdSTEM.qya.^2;
q1 = sqrt(q2);
alpha = q1 * emdSTEM.lambda;
alphaInds = round((alpha + dr/2) / dr);
% alphaInds(alphaInds<1) = 1;
alphaInds(alphaInds>length(emdSTEM.MULTIdetectorAngles)) = 0;
alphaMask = alphaInds > 0;
alphaIndsSub = alphaInds(alphaMask);
Ndet = length(emdSTEM.MULTIdetectorAngles);
emdSTEM.alphaHighAnglesMask = ~alphaMask;

% Initial probe
qProbeMax = emdSTEM.MULTIprobeSemiangleArray / emdSTEM.lambda;
dq = mean([emdSTEM.qxa(2,1) emdSTEM.qya(1,2)]);
PsiProbeInit = (erf((qProbeMax - q1) ...
    /(0.25*dq))*0.5 + 0.5);
PsiProbeInit(:) = PsiProbeInit ...
    .* exp((-1i*pi*emdSTEM.lambda ...
    * emdSTEM.MULTIprobeDefocusArray)*q2);
PsiProbeInit(:) = PsiProbeInit(:) ...
    / sqrt(sum(abs(PsiProbeInit(:)).^2));
emdSTEM.PsiProbeInit = PsiProbeInit;


% Transmission grating
trans = exp(1i*emdSTEM.sigma*emdSTEM.pot(:,:,:,indFP));


% propagators and mask
emdSTEM.qMax = min(max(abs(qx)),max(abs(qy)))/2;
qMask = false(emdSTEM.imageSize);
xOutInds = [(1:(emdSTEM.imageSize(1)/4)) ...
    ((1-emdSTEM.imageSize(1)/4):0)+emdSTEM.imageSize(1)];
yOutInds =  [(1:(emdSTEM.imageSize(2)/4)) ...
    ((1-emdSTEM.imageSize(2)/4):0)+emdSTEM.imageSize(2)];
qMask(xOutInds,yOutInds) = true;
emdSTEM.qMask = qMask;
emdSTEM.prop = qMask ...
    .* exp((-1i*pi*emdSTEM.lambda*emdSTEM.sliceThickness)*q2);


if flagPlot == 1
    xSub = 1:emdSTEM.imageSize(1);
    ySub = 1:emdSTEM.imageSize(2);
end


% init
emdSTEM.MULTIstack = zeros( ...
    length(emdSTEM.MULTIxp),...
    length(emdSTEM.MULTIyp),...
    length(emdSTEM.MULTIdetectorAngles),...
    emdSTEM.numPlanes);
emdSTEM.MULTIhighAngle = zeros( ...
    length(emdSTEM.MULTIxp),...
    length(emdSTEM.MULTIyp),...
    emdSTEM.numPlanes);
if flagKeep4Doutput == true
    emdSTEM.MULTI4D = zeros( ...
        emdSTEM.imageSize(1)/2,...
        emdSTEM.imageSize(2)/2,...
        length(emdSTEM.MULTIxp),...
        length(emdSTEM.MULTIyp));
end

% Main loop
psi = zeros(emdSTEM.imageSize);
% propMask = emdSTEM.prop(emdSTEM.qMask);
% progressbar(0,2);
compCount = compFreq;
emdSTEM.thicknessOutput = (1:emdSTEM.numPlanes)*emdSTEM.sliceThickness;
for a0 = 1:length(emdSTEM.MULTIxp)
    for a1 = 1:length(emdSTEM.MULTIyp)
        % Make probe
        psi(:) = PsiProbeInit ...
            .* exp(-2i*pi ...
            *(emdSTEM.qxa*emdSTEM.MULTIxp(a0) ...
            + emdSTEM.qya*emdSTEM.MULTIyp(a1)));
        
        
        if flagPlot == 1
            Ifft = fftshift(abs(psi));
            Ifft = Ifft(xSub,ySub);
            Ifft = (Ifft - intRangeFFT(1)) ...
                / (intRangeFFT(2) - intRangeFFT(1));
            Ifft = min(max(Ifft,0),1);
            
            I = abs(ifft2(psi));
            I = sqrt(I);
            I = (I(xSub,ySub) - intRange(1)) ...
                / (intRange(2) - intRange(1));
            I = min(max(I,0),1);
            
            if nargin == 2
                fname = [fbase sprintf('%04d',1) '.png'];
                imwrite([I Ifft],fname,'png');
            end
        end
        
        % Propgate through all potential planes
        for a2 = 1:emdSTEM.numPlanes
            %psi = fft2(ifft2(psi).*trans(:,:,a2)).*emdSTEM.prop;
            
            % Keep high angle scattered electrons
            psi = fft2(ifft2(psi).*trans(:,:,a2));
            emdSTEM.MULTIhighAngle(a0,a1,a2) = ...
                emdSTEM.MULTIhighAngle(a0,a1,a2) ...
                + sum(abs(psi(emdSTEM.alphaHighAnglesMask)).^2);
            psi(:) = psi .* emdSTEM.prop;
            
            if flagPlot == 1
                Ifft = fftshift(abs(psi));
                Ifft = Ifft(xSub,ySub);
                Ifft = (Ifft - intRangeFFT(1)) ...
                    / (intRangeFFT(2) - intRangeFFT(1));
                Ifft = min(max(Ifft,0),1);
                
                I = abs(ifft2(psi));
                I = sqrt(I);
                I = (I(xSub,ySub) - intRange(1)) ...
                    / (intRange(2) - intRange(1));
                
                I = min(max(I,0),1);
                
                            if nargin == 2

                fname = [fbase sprintf('%04d',a2+1) '.png'];
                imwrite([I Ifft],fname,'png');
                            end
            end
            
            % save output
            emdSTEM.MULTIstack(a0,a1,:,a2) = ...
                accumarray(alphaIndsSub,abs(psi(alphaMask)).^2,[Ndet 1]);
        end
        if flagKeep4Doutput == true
            emdSTEM.MULTI4D(:,:,a0,a1) = ...
                emdSTEM.MULTI4D(:,:,a0,a1) ...
                + abs(psi(xOutInds,yOutInds)).^2;
        end
        
        
        % Completion
        comp = (a1 / length(emdSTEM.MULTIyp) ...
            + a0 - 1) / length(emdSTEM.MULTIxp);
        %         progressbar(comp,2);
        if comp > compCount
            compOut = round(comp / compFreq) * compFreq;
            % disp([num2str(100*compOut) '% complete'])
            compCount = compCount + compFreq;
        end
    end
end


if flagPlotProbe == true
%     ampRange = [0 5e-6];
    ampRange = [0 2*5e-7];
    ampRangeFFT = [0 1e-2];
    
    pp = exp(1i*angle(psi(1,1)));
    psi = psi / pp;
    
    psiReal = ifft2(psi);
    psiRealplot = colorComplex(abs(psiReal),angle(psiReal),ampRange);
    
    N = size(psi);
    vx = [(1:N(1)/4) ((1-N(1)/4):0)+N(1)];
    vy = [(1:N(2)/4) ((1-N(2)/4):0)+N(2)];
    psiFFT = abs(fftshift(psi(vx,vy)));
    psiFFT(:) = (psiFFT - ampRangeFFT(1)) / (ampRangeFFT(2) - ampRangeFFT(1));
    psiFFT(:) = min(max(psiFFT,0),1);
    psiFFTplot = ind2rgb(round(psiFFT*1023)+1,jetBlack(1024));
    
    figure(3)
    clf
    imagesc(psiRealplot)
    axis equal off
    set(gca,'position',[0 0 1 1])
    
    figure(4)
    clf
    imagesc(psiFFTplot)
    axis equal off
    set(gca,'position',[0 0 1 1])
end


% figure(1)
% clf
% imagesc(fftshift(abs(ifft2(psi)).^2))
% axis equal off
% colormap(gray(256))


% toc
end