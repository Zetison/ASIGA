function data = makeDynamic(data, options, omega)

if ~options.plotTimeOscillation
    return
end
if isfield(options,'noSteps')
    noSteps = options.noSteps;
else
    noSteps = 30;
end
temp = data;
data = zeros(size(temp,1),size(temp,2),noSteps);
for i = 1:noSteps
    t = (i-1)/noSteps*2*pi/omega;
    data(:,:,i) = temp*exp(-1i*omega*t);
end
