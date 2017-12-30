classdef detectorRig < handle
    
    properties
        r;
        N;
        data;
        data_filt;
        data_tof;
        
        data_rectified;
        data_derectified
        
        bp_im;
        tof_im;
        imS;
        
        emissionProbTbl
        absCoeffTbl;
    end
    
    %%
    methods
        %%
        function obj = detectorRig(r,N)
            if r > 0
                obj.r = r;
            end
            if N > 10
                obj.N = N;
            end
            
            obj.data = zeros(obj.N, obj.N);
            obj.data_filt = zeros(obj.N, obj.N);
            
            obj.imS = 100;

            obj.absCoeffTbl = ones(obj.imS) .* 0;
            obj.emissionProb();
            
            obj.data_tof = cell(obj.N);


        end
        
        
        %%
        function detectEmission(obj, loc, theta)
            k = tan(theta);
            b = loc(2) - k*loc(1);

            S2 = [ -( sqrt( (k^2+1)*obj.r^2 - b^2 ) + b*k ) / (k^2 +1),...
                ( sqrt( (k^2+1)*obj.r^2 - b^2 ) - b*k ) / (k^2 +1) ];
            y2 = k.*S2 + b;
           
            angles = angle([S2(1) + y2(1)*1i, S2(2) + y2(2)*1i]);
            
            detectors = floor(angles ./ (2*pi/(obj.N)));
            
            detectors = obj.update_data( detectors );
            if sum(detectors) ~= 0
                obj.detecTimeOfFlight(detectors, loc);
            end
        end

        %%
        function detectors = update_data(obj, detectors)
            if detectors(1) < 0
                detectors(1) = obj.N+detectors(1);
            end
            if detectors(2) < 0
                detectors(2) = obj.N+detectors(2);
            end
            
            if( rand(1) <= obj.emissionProbTbl(...
                                    detectors(1)+1, detectors(2)+1 ) )
                obj.data(detectors(1)+1, detectors(2)+1) = ...
                    obj.data(detectors(1)+1, detectors(2)+1) + 1;
                obj.data(detectors(2)+1, detectors(1)+1) = ...
                    obj.data(detectors(2)+1, detectors(1)+1) + 1;
            else
                detectors = [0,0];
            end
        end
        
        %%
        function filter(obj, data)
            if nargin < 2
                data = obj.data;
            end
            
            fdata = fftshift(fft(data, [], 2), 2);
            h = abs( linspace(-1,1,obj.N) );
                h = repmat(h, [obj.N, 1]);
            
            fdata_filt = ifftshift(fdata .* h, 2);
            obj.data_filt = real( ifft(fdata_filt, [], 2) );

            
            
        end
        
        %% 
        function [pixInd, dist] = pixBetwDet(obj, i, j)
                
                    xi = 1 + round((obj.imS-1) *...
                        (obj.r + ( obj.r *...
                        (cos((i-1)/obj.N * 2*pi)) ))/(2*obj.r));
                    xj = 1 + round((obj.imS-1) *...
                        (obj.r + ( obj.r *...
                        (cos((j-1)/obj.N * 2*pi)) ))/(2*obj.r));
                    
                    yi = 1 + round((obj.imS-1) *...
                        (obj.r + ( obj.r *...
                        (sin((i-1)/obj.N * 2*pi)) ))/(2*obj.r));
                    yj = 1 + round((obj.imS-1) *...
                        (obj.r + ( obj.r *...
                        (sin((j-1)/obj.N * 2*pi)) ))/(2*obj.r));
                    
                    maxDiff = max(abs(xi-xj), abs(yi-yj)) +1 ;
                    
                    pixSub = [round(linspace(xi,xj,maxDiff));...
                        round(linspace(yi,yj,maxDiff))]';
                    pixInd = sub2ind([obj.imS, obj.imS],...
                        pixSub(:,1), pixSub(:,2));
                    
                    dist = sqrt((xi-xj)^2 + (yi-yj)^2);
        end
        
        %%
        function back_project(obj, opt)
            
            if nargin < 2 || ~isfield(opt, 'method')
                opt.method = '';
            end
            
            if strcmp(opt.method, 'unit')
                tau = 1e5;
                ATA = obj.data * obj.data';
                obj.data_filt = pinv(ATA + tau.*eye(size(ATA))) * obj.data;
            elseif strcmp(opt.method, 'none')
                obj.data_filt = obj.data;
            else
                obj.filter(obj.data);
            end
            
            if nargin == 2 && isfield(opt, 'imS')
                obj.imS = opt.imS;
            end
            
            obj.bp_im = zeros(opt.imS);
            for i =  1 : obj.N
                for j = 1 : obj.N-2
                    pixInd = obj.pixBetwDet(i, j);
                    obj.bp_im(pixInd) = obj.bp_im(pixInd ) +...
                        obj.data_filt(i, j);
                end
            end
        end % back_project()
        
        %%
        function rectify_data(obj)
            obj.data_rectified = zeros(obj.N, obj.N);
            for i = 2:obj.N
                obj.data_rectified(i,:) = circshift(obj.data(i,:), -(i-1));
            end
        end
        
        function derectify_data(obj)
            obj.data_derectified = zeros(obj.N, obj.N);
            for i = 2:obj.N
                obj.data_derectified(i,:) = circshift(obj.data_filt(i,:), i-1);
            end
        end
        
        %%
        function emissionProb(obj)
            probTmp = nan(obj.N);
            
            for i = 1 : obj.N-1
                for j = i+1 : obj.N
                    [pixInd, dist] = pixBetwDet(obj, i, j, obj.imS);
                    probTmp(i,j) = exp(-1*...
                        mean(obj.absCoeffTbl(pixInd)) * (dist/obj.imS*2*obj.r));
                end
            end
            
            obj.emissionProbTbl = triu(probTmp) + triu(probTmp)';
            obj.emissionProbTbl( isnan(obj.emissionProbTbl) ) = 0;
            
        end
        

        %%
        function detecTimeOfFlight(obj, detectors, loc)
            c = 2.998e8; % m/s
            
            xd = 1 + round((obj.imS-1) *...
                        (obj.r + ( obj.r *...
                        (cos((detectors-1)/obj.N * 2*pi)) ))/(2*obj.r));
                    
            yd = 1 + round((obj.imS-1) *...
                        (obj.r + ( obj.r *...
                        (sin((detectors-1)/obj.N * 2*pi)) ))/(2*obj.r));
                    
            det_loc = [xd(1), yd(1); xd(2), yd(2)] .* (2*obj.r / obj.imS) - obj.r;
            
            dists = (det_loc - loc).^2;
                dists = sqrt( dists(:,1) + dists(:,2) );
            
            tof = dists ./ (100*c); % s
            
            obj.data_tof{detectors(1)+1, detectors(2)+1} = cat(1,...
                    obj.data_tof{detectors(1)+1, detectors(2)+1},  tof(1));
            obj.data_tof{detectors(2)+1, detectors(1)+1} = cat(1,...
                    obj.data_tof{detectors(2)+1, detectors(1)+1},  tof(2));
        end % detectTimeOfFlight
        
        
        %% 
        function reconstructTimeOfFlight(obj, accuracy)
            obj.tof_im = zeros(obj.imS, obj.imS);
            
            for i = 1 : obj.N-1
                for j = i+1 : obj.N
                    if ~isempty( obj.data_tof(i,j) )
                        
                        pixInd = obj.pixBetwDet(i, j);
                        
                        ratios = obj.data_tof{i,j} ./...
                            (obj.data_tof{i,j} + obj.data_tof{j,i});
                        
                        pixels = pixInd( round(length(pixInd) .* ratios) );
                        
                        obj.tof_im(pixels) = obj.tof_im(pixels) + 1;
                    end                    
                end
            end
            
            
            
        end % reconstructTimeOfFlight
        
        
        
        
        
    end % methods
    
    
    
    
    
end
























