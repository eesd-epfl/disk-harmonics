function writePVD(pvdFilename, base_vtk_files)
    % Function to write ParaView Data (PVD) file for time-dependent data

    % Open the PVD file for writing
    fid = fopen(pvdFilename, 'w');

    % Write PVD header
    fprintf(fid, '<?xml version="1.0"?>\n');
    fprintf(fid, '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid, '  <Collection>\n');

    % Write individual VTK file entries for each time step
    for timestep = 0:numel(base_vtk_files)-1
        vtkFilename = base_vtk_files{timestep+1};
        fprintf(fid, '    <DataSet timestep="%d" group="" part="0" file="%s" name=""/>\n', timestep, vtkFilename);
    end

    % Close the PVD file
    fprintf(fid, '  </Collection>\n');
    fprintf(fid, '</VTKFile>\n');
    fclose(fid);

    disp(['PVD file "', pvdFilename, '" written successfully.']);
end
