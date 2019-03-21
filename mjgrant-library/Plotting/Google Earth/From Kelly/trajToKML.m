function trajToKML(filename, hist)
% filename = 'test.kml';
% 
% hist.longitude = [ -122 -120 -118];
% hist.latitude = [12 14 16];
% hist.altitude = [120000 100000 80000];
 
 
fid = fopen(filename,'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<kml xmlns="http://www.opengis.net/kml/2.2">\n');
fprintf(fid,'<Document>\n');
fprintf(fid,'  <name>Entry</name>\n');
fprintf(fid,'  <open>1</open>\n');
fprintf(fid,'  <Style id="traj">\n');
fprintf(fid,'   <LineStyle>\n');
fprintf(fid,'    <color>7fffff00</color>\n');
fprintf(fid,'    <width>2</width>\n');
fprintf(fid,'    </LineStyle>\n');
fprintf(fid,'   <IconStyle>\n');
fprintf(fid,'    <color>ff00ff7f</color>\n');
fprintf(fid,'    <Icon>\n');
fprintf(fid,'     <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n');
fprintf(fid,'    </Icon>\n');
fprintf(fid,'   </IconStyle>\n');
fprintf(fid,'  </Style>\n');
fprintf(fid,'  <Placemark>\n');
fprintf(fid,'   <styleUrl>#traj</styleUrl>\n');
fprintf(fid,'   <name>Entry Traj</name>\n');
fprintf(fid,'   <LineString>\n');
fprintf(fid,'     <extrude>0</extrude>\n');
fprintf(fid,'     <tessellate>1</tessellate>\n');
fprintf(fid,'     <altitudeMode>relativeToGround</altitudeMode>\n');
fprintf(fid,'     <coordinates>\n');
for i = 1 : 1 : length(hist.latitude)
    fprintf(fid,'%6.4f,%6.4f,%6.4f \n',hist.longitude(i), hist.latitude(i), hist.altitude(i));
end;
fprintf(fid,'     </coordinates>\n');
fprintf(fid,'   </LineString>\n');
fprintf(fid,'  </Placemark>\n');
 
fprintf(fid,'  <Placemark>\n');
fprintf(fid,'   <styleUrl>#traj</styleUrl>\n');
fprintf(fid,'   <name>Entry Traj</name>\n');
fprintf(fid,'   <Point>');
fprintf(fid,'     <coordinates>');
n = length(hist.latitude);
fprintf(fid,'%6.4f,%6.4f,%6.4f \n',hist.longitude(n), hist.latitude(n), hist.altitude(n));
fprintf(fid,'     </coordinates>');
fprintf(fid,'   </Point>');
fprintf(fid,'  </Placemark>\n');
 
fprintf(fid,'</Document>\n');
fprintf(fid,'</kml>\n');
 
fclose(fid);
 
return;
