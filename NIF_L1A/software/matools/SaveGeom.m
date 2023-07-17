UpdateGeom
if nf >0
	save FaultGeom geom -ascii
	disp('New Fault Geometry File Saved')
end

if nv > 0
	save VolGeom vgeom -ascii
	disp('New Volume Geometry File Saved')
end
