function writeVTU(vertices, faces, dataStructArray, filename)
% writeVTU - Function to write VTU file in ASCII format with multiple data structures
%
%   writeVTU(vertices, faces, dataStructArray, filename)
%
%   Parameters:
%       - vertices: Vertex coordinates, N x 3 matrix
%       - faces: Face indices, M x 3 matrix (triangular mesh)
%       - dataStructArray: Array of structures containing additional data
%         Each structure has fields: name (string), data (N x D matrix)
%       - filename: Output VTU file name (including path if needed)
%
%   Example:
%       vertices = [0 0 0; 1 0 0; 0 1 0; 1 1 0];
%       faces = [1 2 3; 2 3 4];
%       dataStructArray(1).name = 'ScalarData';
%       dataStructArray(1).data = [0.1; 0.2; 0.3; 0.4];
%       writeVTU(vertices, faces, dataStructArray, 'output.vtu');
%
%   Credits:
%       This function is based on collaborative work with OpenAI's GPT-3 model.
%       Assistance provided by ChatGPT, a product of OpenAI.
%
%   Note: Mahmoud Shaqfa modified this file that was primitively provided by
%   ChatGPT-3 @ 2024.

% Open the file for writing
fid = fopen(filename, 'w');

% Write VTU header
fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid, '  <UnstructuredGrid>\n');
fprintf(fid, '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n', size(vertices, 1), size(faces, 1));

% Write vertex information
fprintf(fid, '      <Points>\n');
fprintf(fid, '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
fprintf(fid, '          %f %f %f\n', vertices');
fprintf(fid, '        </DataArray>\n');
fprintf(fid, '      </Points>\n');

% Write cell information
fprintf(fid, '      <Cells>\n');
fprintf(fid, '        <DataArray type="Int32" Name="connectivity" format="ascii">\n');
fprintf(fid, '          %d\n', reshape(faces' - 1, 1, []));
fprintf(fid, '        </DataArray>\n');
fprintf(fid, '        <DataArray type="Int32" Name="offsets" format="ascii">\n');
fprintf(fid, '          %d\n', 3 * (1:size(faces, 1)));
fprintf(fid, '        </DataArray>\n');
fprintf(fid, '        <DataArray type="UInt8" Name="types" format="ascii">\n');
fprintf(fid, '          %d\n', repmat(5, size(faces, 1), 1));  % 5 corresponds to VTK_TRIANGLE
fprintf(fid, '        </DataArray>\n');
fprintf(fid, '      </Cells>\n');

% Write data
fprintf(fid, '      <PointData>\n');

for i = 1:numel(dataStructArray)
    dataStruct = dataStructArray(i);
    dataName = dataStruct.name;
    data = dataStruct.data;
    dataDim = size(data, 2);

    fprintf(fid, '        <DataArray type="Float32" Name="%s" NumberOfComponents="%d" format="ascii">\n', dataName, dataDim);
    fprintf(fid, '          %f\n', data');
    fprintf(fid, '        </DataArray>\n');

end

fprintf(fid, '      </PointData>\n');

% Close the VTU file
fprintf(fid, '    </Piece>\n');
fprintf(fid, '  </UnstructuredGrid>\n');
fprintf(fid, '</VTKFile>\n');

fclose(fid);
end

