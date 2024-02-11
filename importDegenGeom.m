function panXY = importDegenGeom( dgfile )
%importDgenGeom Reads OpenVSP degengeom file
%   Detailed explanation goes here

  panXY = vspPlotDegenStick(dgfile);
  
  vspPlotDegenPlate(dgfile);
  
  vspPlotDegenSurf(dgfile);
end

