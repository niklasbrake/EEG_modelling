classdef clrsPT
    % Paul tol vibrant colour scehem
    properties (Constant)
        qualitative_CM = struct('blue',[0,119,187]/255, ...
                            'cyan',[51,187,238]/255, ...
                            'teal',[0,153,136]/255, ...
                            'orange',[238,119,51]/255, ...
                            'red',[204,51,17]/255, ...
                            'magenta',[238,51,119]/255, ...
                            'grey',[187,187,187]/255);
        diverging_CM = [[33,102,172]; ...
                    [67,147,195]; ...
                    [146,197,222]; ...
                    [209,229,240]; ...
                    [247,247,247]; ...
                    [253,219,199]; ...
                    [244,165,130]; ...
                    [214,96,77]; ...
                    [178,24,43]]/255;
        iridescent_CM = [254,251,233; ...
                      252,247,213; ...
                      245,243,193;
                      234,240,181;
                      221,236,191;
                      208,231,202;
                      194,227,210;
                      181,221,216;
                      168,216,220;
                      155,210,225;
                      141,203,228;
                      129,196,231;
                      123,188,231;
                      126,178,228;
                      136,165,221;
                      147,152,210;
                      155,138,196;
                      157,125,178;
                      154,112,158;
                      144,99,136;
                      128,87,112;
                      104,73,87;
                      70,53,58]/255;
            sequential_CM = [255,255,255; ...
                    255,247,188; ...
                    254,227,145; ...
                    254,196,79; ...
                    251,154,41; ...
                    236,112,20; ...
                    204,76,2; ...
                    153,52,3; ...
                    102,37,6]/255;
    end
    methods (Static)
        function clrs = lines(N)
            colourOrder = [1,4,2,5,3,6,7];
            temp = clrsPT.qualitative_CM;
            fn = fieldnames(temp);
            for i = 1:length(fn)
                clrs(i,:) = temp.(fn{i});
            end
            if(nargin==0)
                clrs = clrs(colourOrder,:);
            else
                colourOrder = repmat(colourOrder,1,ceil(N/length(colourOrder)));
                clrs = clrs(colourOrder(1:N),:);
            end
        end
        function clrs = diverging(N)
            I1 = linspace(0,1,9);
            I2 = linspace(0,1,N);
            clrs = interp1(I1,clrsPT.diverging_CM,I2);
        end
        function clrs = sequential(N)
            I1 = linspace(0,1,9);
            I2 = linspace(0,1,N);
            clrs = interp1(I1,clrsPT.sequential_CM,I2);
        end
        function clrs = iridescent(N)
            I1 = linspace(0,1,23);
            I2 = linspace(0,1,N);
            clrs = interp1(I1,clrsPT.iridescent_CM,I2);
        end
    end
end
