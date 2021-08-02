function structDP = getImportedDP(file)

% Clarissa Towle - 2020 June

% This script obtains the diffraction pattern for an Au nanoparticle or an 
%   MoS2 film of a given layer thickness.

% Note all variables are in Angstroms for distances, rads for angles

% Inputs:
% file -- the .csv file of atom positions and identities

% Outputs:
% structDP -- simulated diffraction pattern based on atoms in 'file'

% Input variables:
atoms = readmatrix(file);
flagPlot = true;
cellDim = [250 250 30];


% Plotting to show output
if flagPlot == true
    figure(11)
    clf
    set(gcf,'color','w')
    hold on
    s = atoms(:,4) == 12;
    scatter3(atoms(s,2),atoms(s,1),atoms(s,3),...
        'marker','o','sizedata',8,'linewidth',0.5,...
        'markeredgecolor','none','markerfacecolor',[1 1 1])
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
    
%     % UC boundaries
%     line([0 0 cellDim(2) cellDim(2) 0],...
%         [0 cellDim(1) cellDim(1) 0 0],...
%         [0 0 0 0 0],'linewidth',2,'color','k')
%     line([0 0 cellDim(2) cellDim(2) 0],...
%         [0 cellDim(1) cellDim(1) 0 0],...
%         [0 0 0 0 0]+cellDim(3),'linewidth',2,'color','k')
    
    hold off
    
    axis equal on
%     view ([0 0 1])
    view([-5 -2 5])
    set(gca,'BoxStyle','full','Box','on')
%     set(gca,'position',[0 0 1 1],'ydir','reverse')
    xlim([0 cellDim(2)]+[-1 1])
    ylim([0 cellDim(1)]+[-1 1])
%     zlim([0 cellDim(3)]+[-1 1])
    
end

emdSTEM = PRISM01(atoms,cellDim);

emdSTEM = PRISMmultislice(emdSTEM);

% plotImage(fftshift(emdSTEM.MULTI4D),'absolute',jetBlack,[0 0.00005]);
structDP=fftshift(emdSTEM.MULTI4D);

end
