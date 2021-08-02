function avgInt = getMeanIntensity2(species,nlayers)

% Colin Ophus - 2020 Feb
% Modified Clarissa Towle - 2020 April

% This script generates the shape and structure of a single nanoparticle,
% rotates it to a given orientation, and then plonks it into a substrate.

% This script has been modified to obtain the mean intensity for a set of
% diffraction disks belonging to either Au or MoS2, in a sample that is an
% Au nanoparticle or an MoS2 film of a given layer thickness.

% Run over an array of diffraction patterns (4D STEM), this script produces
% a type of dark-field image.

% Note all variables are in Angstroms for distances, rads for angles

% Outputs:
% atoms - [N x 4] array where columns correspond to [x y z atomic_number]
% cellDim - 3 element vector [x_cell_dim y_cell_dim z_cell_dim]

% Input variables:
flagPlot = true;
substrateSizeApprox = [1 1]*150;
substrateNumLayers = (nlayers-1)/2;

shiftSubstrateZ = 2.3;  % distance from base of NP to substrate
shiftNP = [0 -0.83];  % x,y shifting of the NP
dzBound = 1;  % Extra spacing at top and bottom of sim cell

% Facet type direction, then length in Angstroms
shapeNP = [ ...
    1 1 1 50;
    1 0 0 100];
% shapeNP = [ ...
%     1 1 1 38;
%     1 0 0 80];

% Distance from NP origin to substrate.
% Controls how much NP wets the surface.
if strcmpi(species,'au')
    distSubstrateNP = -49+2.3545*(nlayers-1);
%     distSubstrateNP = -37+2.3545*(nlayers-1);
elseif strcmpi(species,'mos2')
    distSubstrateNP = -49;
%     distSubstrateNP = -37;
end


% Euler rotation angles around z x z axis in this order
zxzRotationNP = [45 54.7356 0] * pi / 180;% 111 orientation
% zxzRotationNP = [0 45 90] * pi / 180;  % 110 orientation
% zxzRotationNP = [0 0 45] * pi / 180;  % 100 orientation

% zxzRotationNP = [0 71.5651 90] * pi / 180;  % 310 orientation


% Structural data:
% Lattice parameter and atom basis for gold:
aNP = 4.0782;
basisNP = [ ...
    0.0 0.0 0.0 79;
    0.5 0.5 0.0 79;
    0.5 0.0 0.5 79;
    0.0 0.5 0.5 79];
% Lattice parameter and atom basis for top layer MoS2
abcSub1 = [2.8837 2.8837*sqrt(3) 12.3];
d = 0.121;
basisSub1 = [ ...
    0/2 1/6 1/4 42;
    1/2 4/6 1/4 42;
    0/2 5/6 3/4 42;
    1/2 2/6 3/4 42;
    ...
    0/2 1/6 3/4+d 16;
    0/2 1/6 3/4-d 16;
    1/2 4/6 3/4+d 16;
    1/2 4/6 3/4-d 16;
    0/2 5/6 1/4+d 16;
    0/2 5/6 1/4-d 16;
    1/2 2/6 1/4+d 16;
    1/2 2/6 1/4-d 16;
    ];
% Lattice parameter and atom basis for rest of MoS2
abcSub = [3.1604 3.1604*sqrt(3) 12.3];
d = 0.121;
basisSub = [ ...
    0/2 1/6 1/4 42;
    1/2 4/6 1/4 42;
    0/2 5/6 3/4 42;
    1/2 2/6 3/4 42;
    ...
    0/2 1/6 3/4+d 16;
    0/2 1/6 3/4-d 16;
    1/2 4/6 3/4+d 16;
    1/2 4/6 3/4-d 16;
    0/2 5/6 1/4+d 16;
    0/2 5/6 1/4-d 16;
    1/2 2/6 1/4+d 16;
    1/2 2/6 1/4-d 16;
    ];



% Generate substrate
numXYsubstrate = round(substrateSizeApprox ./ abcSub(1:2));
cellDimXY = numXYsubstrate .* abcSub(1:2);

disp(size(cellDimXY))

% Tile out unit cells
[ya,xa,za] = meshgrid( ...
    (0:numXYsubstrate(2)-1),(0:numXYsubstrate(1)-1),0:ceil(substrateNumLayers-1));
p = [xa(:) ya(:) za(:)];
[bInd1,pInd1] = meshgrid(1:size(basisSub1,1),1:size(p,1));
[bInd,pInd] = meshgrid(1:size(basisSub,1),1:size(p,1));
atomsSub1 = [p(pInd1(:),:) + basisSub1(bInd1(:),1:3) basisSub1(bInd1(:),4)];
atomsSub = [p(pInd(:),:) + basisSub(bInd(:),1:3) basisSub(bInd(:),4)];
% Remove extra half plane if needed
del1 = atomsSub1(:,3) > 0.5;
atomsSub1(del1,:) = [];
del = atomsSub(:,3) > substrateNumLayers;
atomsSub(del,:) = [];
% Scale to Angstroms
atomsSub1(:,1:3) = atomsSub1(:,1:3) .* abcSub1;
atomsSub(:,1:3) = atomsSub(:,1:3) .* abcSub;



% Generate nanoparticle block
maxRadius = round(max(sqrt(3) * shapeNP(:,4)) / aNP);
v = -maxRadius:maxRadius;
[ya,xa,za] = meshgrid(v,v,v);
p = [xa(:) ya(:) za(:)];
[bInd,pInd] = meshgrid(1:size(basisNP,1),1:size(p,1));
atomsNP = [p(pInd(:),:) + basisNP(bInd(:),1:3) basisNP(bInd(:),4)];
% Scale to Angstroms
atomsNP(:,1:3) = atomsNP(:,1:3) * aNP;



% Remove all Wulff facets from NP
for a0 = 1:size(shapeNP,1)
    % Facet data
    n = shapeNP(a0,1:3);
    n = n / norm(n);
    dist = shapeNP(a0,4);
    
    % Perms to generate all permutations
    %     nAll = perms(perms(n,'signs'),'unique');
    nAll = permsColin(n,'signs','unique');
    
    % Cut facets from NP block
    for a1 = 1:size(nAll,1)
        del = atomsNP(:,1)*nAll(a1,1) ...
            + atomsNP(:,2)*nAll(a1,2) ...
            + atomsNP(:,3)*nAll(a1,3) ...
            > dist;
        atomsNP(del,:) = [];
    end
end

% NP rotation
% 1st rotation - z axis
theta = zxzRotationNP(1);
m = [cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];
atomsNP(:,1:3) = atomsNP(:,1:3) * m;
% 2nd rotation - x axis
theta = zxzRotationNP(2);
m = [1 0 0;
    0 cos(theta) -sin(theta);
    0 sin(theta) cos(theta)];
atomsNP(:,1:3) = atomsNP(:,1:3) * m;
% 3rd rotation - z axis
theta = zxzRotationNP(3);
m = [cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];
atomsNP(:,1:3) = atomsNP(:,1:3) * m;


% Substrate wetting / deletion
del = atomsNP(:,3) < -distSubstrateNP;
atomsNP(del,:) = [];


% Top layer MoS2 rotation - z axis
theta1 = 180*pi/180;
m1 = [cos(theta1) -sin(theta1) 0;
    sin(theta1) cos(theta1) 0;
    0 0 1];
atomsSub1(:,1:3) = atomsSub1(:,1:3) * m1;


% Assemble entire unit cell
atomsSub1(:,3) = atomsSub1(:,3) - min(atomsSub1(:,3));
atomsSub(:,3) = atomsSub(:,3) - min(atomsSub(:,3));
atomsNP(:,3) = atomsNP(:,3) - min(atomsNP(:,3));
zSub = 2*max(atomsSub1(:,3))+max(atomsSub(:,3));
zNP = max(atomsNP(:,3));
% Shift atoms along beam direction
atomsSub(:,3) = atomsSub(:,3) + dzBound;
atomsSub1(:,3) = atomsSub1(:,3) + dzBound + zSub - max(atomsSub1(:,3));

disp(size(atomsNP(:,3)))
disp(size(dzBound))
disp(size(zSub))
disp(size(shiftSubstrateZ))

% zSub
atomsNP(:,3) = atomsNP(:,3) + dzBound + shiftSubstrateZ;
% Move MoS2 top layer to center of cell
atomsSub1(:,1) = atomsSub1(:,1) + 0.95*cellDimXY(1);
atomsSub1(:,2) = atomsSub1(:,2) + 0.95*cellDimXY(2);
% Move NP to center of cell
atomsNP(:,1) = atomsNP(:,1) + shiftNP(1) + cellDimXY(1)/2;
atomsNP(:,2) = atomsNP(:,2) + shiftNP(2) + cellDimXY(2)/2;


% put everything together
if nlayers > 1
    if strcmpi(species, 'mos2')
        atoms = [atomsSub1; atomsSub];
    elseif strcmpi(species, 'both')
        atoms = [atomsSub1; atomsSub; atomsNP];
    else atoms = [atomsSub1; atomsSub; atomsNP];
    end
elseif nlayers == 1
    if strcmpi(species, 'mos2')
        atoms = [atomsSub1];
    elseif strcmpi(species, 'both')
        atoms = [atomsSub1; atomsNP];
    else atoms = [atomsSub1; atomsNP];
    end
else
    disp('Error: There has to be at least 1 layer. Please enter a positive integer.')
end

cellDimZ = 2*dzBound + shiftSubstrateZ + zSub + zNP;
cellDim = [cellDimXY cellDimZ];

% Plotting to show output
if flagPlot == true
    figure(11)
    clf
    set(gcf,'color','w')
    hold on
    s = atoms(:,4) == 16;
    scatter3(atoms(s,2),atoms(s,1),atoms(s,3),...
        'marker','o','sizedata',10,'linewidth',0.5,...
        'markeredgecolor','none','markerfacecolor',[1 0.576 0.576])
    s = atoms(:,4) == 42;
    scatter3(atoms(s,2),atoms(s,1),atoms(s,3),...
        'marker','o','sizedata',15,'linewidth',0.5,...
        'markeredgecolor',[0 0 0],'markerfacecolor',[0.333 0.498 1])
    s = atoms(:,4) == 79;
    scatter3(atoms(s,2),atoms(s,1),atoms(s,3),...
        'marker','o','sizedata',20,'linewidth',0.5,...
        'markeredgecolor',[0 0 0],'markerfacecolor',[1 0.92 0.153])
    
    disp(size(cellDim))
    
    % UC boundaries
    line([0 0 cellDim(2) cellDim(2) 0],...
        [0 cellDim(1) cellDim(1) 0 0],...
        [0 0 0 0 0],'linewidth',2,'color','k')
    line([0 0 cellDim(2) cellDim(2) 0],...
        [0 cellDim(1) cellDim(1) 0 0],...
        [0 0 0 0 0]+cellDim(3),'linewidth',2,'color','k')
    
    hold off
    
    axis equal off
    view([0 0 1])
    set(gca,'position',[0 0 1 1],'ydir','reverse')
    xlim([0 cellDim(2)]+[-1 1])
    ylim([0 cellDim(1)]+[-1 1])
    zlim([0 cellDim(3)]+[-1 1])
    
end

emdSTEM = PRISM01(atoms,cellDim);

emdSTEM = PRISMmultislice(emdSTEM);

plotImage(fftshift(emdSTEM.MULTI4D),'absolute',jetBlack,[0 0.0002]);
structDP=fftshift(emdSTEM.MULTI4D);


% For this sample, the coordinate ranges for a set of Au diffraction disks
% and a set of MoS2 diffraction disks is hard-coded below. Replace with
% the specified disks desired.
auspots={312:323 334:345; 312:323 393:404; 415:426 334:345; 415:426 393:404; 364:374 305:315; 364:374 423:433};
% mos2spots={318:326 338:346; 318:326 392:400; 412:420 392:400; 412:420 338:346; 365:373 311:319; 365:373 419:427};

y = [];
for t = 1:6

% if strcmpi(species,'mos2')
%     x = sum(structDP(mos2spots{t,1:2}),'all');
% elseif strcmpi(species,'au')
x = sum(structDP(auspots{t,1:2}),'all');
% end

y = [y ; x];
end
avgInt = mean(y(1:6),'all');

end