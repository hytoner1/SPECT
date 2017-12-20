

dR = detectorRig(5, 200);
locs = rand(1e5,2) .* 2 - 1;
phis = rand(1e5,1) .* pi;
% locs = zeros(1000, 2);
% phis = linspace(0,2*pi,1000);

figure; 
    scatter(locs(:,1), locs(:,2),'.');
    viscircles([0,0], dR.r);
   
%%
    
for i=1:size(locs,1)
    dR.detectEmission(locs(i,:), phis(i))
end

figure;
    imagesc(dR.data)
    
    
%%
opt.method = 'unit';
opt.imS = 50;
dR.back_project(opt);

figure;
    imagesc(dR.bp_im );
    ax = gca;
    ax.XTick = linspace(0, size(dR.bp_im,1), 11);
        ax.XTickLabel = linspace(-dR.r, dR.r, 11);
    ax.YTick = linspace(0, size(dR.bp_im,1), 11);
        ax.YTickLabel = linspace(-dR.r, dR.r, 11);
        
    %rectangle('Position', [30,30,40,40])
    
