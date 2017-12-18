classdef detectorRig < handle
    
    properties
        r;
        N;
        data;
        data_filt;
        bp_im;
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
            
            obj.initData
        end
        
        %%
        function initData(obj)
            obj.data = zeros(obj.N, obj.N);
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
            
            obj.update_data( detectors );
        end

        %%
        function update_data(obj, detectors)
            if detectors(1) < 0
                detectors(1) = obj.N+detectors(1);
            end
            if detectors(2) < 0
                detectors(2) = obj.N+detectors(2);
            end
            if detectors(1) > detectors(2)
                tmp = detectors(2);
                detectors(2) = detectors(1);
                detectors(1) = tmp;
            end
            
            obj.data(detectors(1)+1, detectors(2)+1) = ...
                obj.data(detectors(1)+1, detectors(2)+1) + 1;
        end
        
        %%
        function back_project(obj, data)
            if nargin < 2
                data = obj.data_filt;
            end
            
            imS = 100;
            obj.bp_im = zeros(imS);
                        
            for i =  1 : obj.N-20
                for j = i+1 : obj.N
                    if(data(i,j) > 0)
                        xi = 1 + round((imS-1) *...
                            (obj.r + ( obj.r *...
                            (cos((i-1)/obj.N * 2*pi)) ))/(2*obj.r));
                        xj = 1 + round((imS-1) *...
                            (obj.r + ( obj.r *...
                            (cos((j-1)/obj.N * 2*pi)) ))/(2*obj.r));
                        
                        yi = 1 + round((imS-1) *...
                            (obj.r + ( obj.r *...
                            (sin((i-1)/obj.N * 2*pi)) ))/(2*obj.r));
                        yj = 1 + round((imS-1) *...
                            (obj.r + ( obj.r *...
                            (sin((j-1)/obj.N * 2*pi)) ))/(2*obj.r));
                        
                        maxDiff = max(abs(xi-xj), abs(yi-yj)) +1 ;
                        
                        asd = [round(linspace(xi,xj,maxDiff));...
                            round(linspace(yi,yj,maxDiff))]';
                        asd = sub2ind([imS, imS], asd(:,1), asd(:,2));
                        
                        obj.bp_im(asd) = obj.bp_im(asd) + data(i,j);
                        
                    end
                end
            end
        end % back_project()
        
        %%
        function filter(obj)
            fdata = fftshift(fft(obj.data, [], 1), 1);
            size(fdata)
            h = abs( linspace(-1,1,obj.N) );
            size(h)
                h = repmat(h, [obj.N, 1]);
            size(h)
            
            fdata_filt = ifftshift( fdata .* h, 1 );
            obj.data_filt = real( ifft(fdata_filt, [], 1) );
            
        end
        
        
    end % methods
    
    
    
    
    
end