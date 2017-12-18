

dR = detectorRig(5, 200);
locs = rand(1000,2) .* 4 - 2;
phis = rand(1000,1) .* pi;

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

dR.back_project();

%%

figure;
    imagesc(sqrt(dR.bp_im));
