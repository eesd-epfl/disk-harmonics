function writeVTK(vertices, faces, dataStructArray, filename)
    % Function to write VTK file in ASCII format with multiple data structures

    % Open the file for writing
    fid = fopen(filename, 'w');

    % Write VTK header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'VTK file written from MATLAB\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET POLYDATA\n');

    % Write vertex information
    numVertices = size(vertices, 1);
    fprintf(fid, 'POINTS %d float\n', numVertices);
    fprintf(fid, '%f %f %f\n', vertices');

    % Write face information
    numFaces = size(faces, 1);
    numFaceVertices = size(faces, 2);
    fprintf(fid, 'POLYGONS %d %d\n', numFaces, numFaces * (numFaceVertices + 1));
    fprintf(fid, '%d %d %d %d\n', [numFaceVertices * ones(numFaces, 1), faces-1]');

    % Write data
    fprintf(fid, 'POINT_DATA %d\n', numVertices);

    for i = 1:numel(dataStructArray)
        dataStruct = dataStructArray(i);
        dataName = dataStruct.name;
        data = dataStruct.data;
        dataDim = size(data, 2);

        if dataDim == 1
            fprintf(fid, 'SCALARS %s float %d\n', dataName, dataDim);
            fprintf(fid, 'LOOKUP_TABLE default\n');
            fprintf(fid, '%f\n', data');
        elseif dataDim == 3
            fprintf(fid, 'VECTORS %s float\n', dataName);
            fprintf(fid, '%f %f %f\n', data');
        else
            error('Unsupported data dimension. Only scalar or 3D vector data is supported.');
        end
    end

    % Close the file
    fclose(fid);

    disp(['VTK file "', filename, '" written successfully.']);
end

% Example usage for multiple data structures
% vertices = [0, 0, 0; 1, 0, 0; 0.5, 1, 0];
% faces = [1, 2, 3];
% scalarData = [1.0; 2.0; 3.0];
% vectorData = [0.1, 0.2, 0.3; 0.4, 0.5, 0.6; 0.7, 0.8, 0.9];
% 
% dataStructArray = [
%     struct('name', 'ScalarData', 'data', scalarData);
%     struct('name', 'VectorData', 'data', vectorData)
% ];
% 
% writeVTK(vertices, faces, dataStructArray, 'output_multiple.vtk');