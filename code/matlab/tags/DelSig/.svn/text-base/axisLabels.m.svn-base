function s = axisLabels(range,incr)
%function s = axisLabels(range,incr)
range(abs(range)<1e-6) = 0;
s = cell(1,length(range));
for i=1:incr:length(range)
    s{i}=sprintf('%g',range(i));
end
