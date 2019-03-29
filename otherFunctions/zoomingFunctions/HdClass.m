classdef HdClass < handle
   properties
      data
   end
   methods
      function obj = HdClass(val)
         if nargin > 0
            obj.data = val;
         end
      end
   end
end