classdef objMatrix   
    properties
        value
    end   
    methods        
        function obj = objMatrix(n,m)
            if nargin > 0
                obj = repmat(obj,n,m);
            end
        end        
    end
end