

dR = detectorRig(5, 200);
% locs = rand(1e4,2) .* [2,1]- 0.5;
% locs = [locs; rand(1e4,2) .* [1, 2] - 0.5];
% phis = rand(size(locs,1),1) .* pi;
% locs = zeros(1000, 2);
% phis = linspace(0,2*pi,1000);

% figure; 
%     scatter(locs(:,1), locs(:,2),'.');
%     viscircles([0,0], dR.r);
   
%%
    
for i=1:1e5 %size(locs,1)
    loc = rand(1,2) .* [2,1]- 0.5;
    loc = [loc; rand(1,2) .* [1, 2] - 0.5];
    phi = rand(1) .* pi;
    
    if rand(1) > 0.5
        loc = loc(1,:);
    else
        loc = loc(2,:);
    end
    dR.detectEmission(loc, phi);
%     dR.detectEmission(locs(i,:), phis(i))
    
end

figure;
    imagesc(dR.data)
    
    
%%
opt.method = '';
opt.imS = 100;
dR.back_project(opt);

figure;
    imagesc(dR.bp_im );
    ax = gca;
    ax.XTick = linspace(0, size(dR.bp_im,1), 11);
        ax.XTickLabel = linspace(-dR.r, dR.r, 11);
    ax.YTick = linspace(0, size(dR.bp_im,1), 11);
        ax.YTickLabel = linspace(-dR.r, dR.r, 11);
        
    %rectangle('Position', [30,30,40,40])
    
%%

dR.reconstructTimeOfFlight();

figure; 
    imagesc( dR.tof_im );
    ax = gca;
    ax.XTick = linspace(0, size(dR.tof_im,1), 11);
        ax.XTickLabel = linspace(-dR.r, dR.r, 11);
    ax.YTick = linspace(0, size(dR.tof_im,1), 11);
        ax.YTickLabel = linspace(-dR.r, dR.r, 11);
        
        
        
        
%% test

im = zeros(100);
figure;
for i=1:200
    im = zeros(100);
    for j=1:200
        pI = dR.pixBetwDet(i,j);
        im(pI) = 1;
        
    end
    imagesc(im);
    pause(0.2);
end

