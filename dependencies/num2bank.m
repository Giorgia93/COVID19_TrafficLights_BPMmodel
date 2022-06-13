function [str] = num2bank(num)
   str = arrayfun(@(x) num2bankScalar(x), num, 'UniformOutput', false) ;
end