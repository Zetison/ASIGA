function makeVTKfile_new_mesh(data, newOptions)


options = struct('name','untitled',...
                 'plotDisplacementVectors',0,...
                 'plotTimeOscillation', 0);

%# read the acceptable names
optionNames = fieldnames(options);

%# count arguments
nArgs = length(newOptions);
if round(nArgs/2) ~= nArgs/2
   error('Must have propertyName/propertyValue pairs')
end

for pair = reshape(newOptions,2,[]) %# pair is {propName;propValue}
   inpName = pair{1}; %# make case insensitive

   if any(strcmp(inpName,optionNames))
      %# overwrite options. If you want you can test for the right class here
      %# Also, if you find out that there is an option you keep getting wrong,
      %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end




vtfFileName = options.name;
nodes = data.nodes;
visElements = data.visElements;
if options.plotTimeOscillation
    noSteps = 30;
    omega = data.omega;
else
    noSteps = 1;
    omega = 1;
end

%--------------------------------------------------------------------------
% Write the header for a structured grid:
%--------------------------------------------------------------------------
noVisElems = length(visElements);
noNodes = size(nodes,1);


parfor (s = 1:noSteps, 321)
% for s = 1:noSteps
    t = (s-1)/noSteps*2*pi/omega;


    if noSteps > 1
        fid = fopen([vtfFileName 'time_' num2str(s) '.vtk'],'wt+','b');
    else
        fid = fopen([vtfFileName '.vtk'],'wt+','b');
    end
    
    fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid,'\t<UnstructuredGrid>\n');
    fprintf(fid,'\t\t<Piece NumberOfPoints="%d" NumberOfCells="%d">\n',noNodes,noVisElems);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Add coordinates
    fprintf(fid,'\t\t\t<Points>\n');
    fprintf(fid,'\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
    
    for i = 1:noNodes
        fprintf(fid,'\t\t\t\t\t%.15g %.15g %.15g\n', [nodes(i,1) nodes(i,2) nodes(i,3)]);
    end
    fprintf(fid,'\t\t\t\t</DataArray>\n');
    fprintf(fid,'\t\t\t</Points>\n');


    fprintf(fid,'\t\t\t<Cells>\n');
    fprintf(fid,'\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">\n');
    for eVis = 1:noVisElems
        fprintf(fid,'\t\t\t\t\t');
        for i = 1:length(visElements{eVis})
            fprintf(fid,'%d ', visElements{eVis}(i)-1);
        end
        fprintf(fid,'\n');
    end
            
    fprintf(fid,'\t\t\t\t</DataArray>\n');

    fprintf(fid,'\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">\n');
    fprintf(fid,'\t\t\t\t\t');
    count = 0;
    for eVis = 1:noVisElems
        count = count + length(visElements{eVis});
        fprintf(fid,'%d ', count);
    end
    fprintf(fid,'\n');
    fprintf(fid,'\t\t\t\t</DataArray>\n');

    fprintf(fid,'\t\t\t\t<DataArray type="UInt8" Name="types" format="ascii">\n');
    fprintf(fid,'\t\t\t\t\t');
    for eVis = 1:noVisElems
        fprintf(fid,'4 ');
    end
    fprintf(fid,'\n');
    fprintf(fid,'\t\t\t\t</DataArray>\n');
    
    fprintf(fid,'\t\t\t</Cells>\n');

    %% Add Displacement (vectors)
    if options.plotDisplacementVectors   
        fprintf(fid,'\t\t\t<PointData Scalars="scalars">\n'); 
        displacement = data.displacement;
        fprintf(fid,'\t\t\t\t<DataArray type="Float32" Name="Displacement" NumberOfComponents="3" format="ascii">\n');
        for i = 1:noNodes
            fprintf(fid,'\t\t\t\t\t%.15g %.15g %.15g\n', real(displacement(i,:)*exp(-1i*omega*t)));
        end
        fprintf(fid,'\t\t\t\t</DataArray>\n');
        fprintf(fid,'\t\t\t</PointData>\n');
    end

    fprintf(fid,'\t\t</Piece>\n');
    fprintf(fid,'\t</UnstructuredGrid>\n');
    fprintf(fid,'</VTKFile>\n');        
    
    fclose(fid);
end
