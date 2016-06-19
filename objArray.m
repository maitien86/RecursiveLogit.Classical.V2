classdef objArray   
    properties
        value
    end   
    methods        
        function obj = objArray(n)
            if nargin > 0
                obj = repmat(obj,1,n);
            end
        end        
    end
end