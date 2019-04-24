
x = nodes(:,1);
y = nodes(:,2);
z = nodes(:,3);

scalarField = u_analytic(x,y,z);


dataAnal.nodes = nodes;
dataAnal.visElements = visElements;
if takeRealPart
    dataAnal.scalarField = real(scalarField);
else
    dataAnal.scalarField = imag(scalarField);
end


options = {'name',vtfFileName,...
     'plotScalarField',1};

makeVTFfile_new(dataAnal, options);
