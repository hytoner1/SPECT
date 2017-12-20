

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
dR.filter();

figure;
    imagesc(dR.data_filt)

%%
% TODO: Add option to use either (f)BP, unit matrix inv or lapl. inv
opt.method = 'unit';
dR.back_project();

figure;
    imagesc(dR.bp_im );
    ax = gca;
    ax.XTick = 0:10:100;
        ax.XTickLabel = -5 : 1 : 5;
    ax.YTick = 0:10:100;
        ax.YTickLabel = -5 : 1 : 5;
        
    %rectangle('Position', [30,30,40,40])
    
