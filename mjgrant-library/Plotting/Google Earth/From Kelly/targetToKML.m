function targetToKML(filename, target)
 
fid = fopen(filename,'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<kml xmlns="http://www.opengis.net/kml/2.2">\n');
fprintf(fid,'<Document>\n');
fprintf(fid,'  <name>Targeted Landing Site</name>\n');
fprintf(fid,'  <open>1</open>\n');
fprintf(fid,'  <Style id="targ">\n');
fprintf(fid,'   <LabelStyle>\n');
fprintf(fid,'    <color>ff00ffff</color>\n');
fprintf(fid,'   </LabelStyle>\n');
fprintf(fid,'   <IconStyle>\n');
fprintf(fid,'    <color>ff00ffff</color>\n');
fprintf(fid,'    <Icon>');
fprintf(fid,'     <href>http://maps.google.com/mapfiles/kml/shapes/target.png</href>');
fprintf(fid,'    </Icon>');
fprintf(fid,'   </IconStyle>\n');
fprintf(fid,'  </Style>\n');
fprintf(fid,'  <Placemark>\n');
fprintf(fid,'   <styleUrl>#targ</styleUrl>\n');
fprintf(fid,'   <name>Target</name>\n');
fprintf(fid,'   <Point>\n');
fprintf(fid,'     <altitudeMode>relativeToGround</altitudeMode>\n');
fprintf(fid,'     <coordinates>\n');
fprintf(fid,'%6.4f,%6.4f,%6.4f \n',target.longitude, target.latitude, 0);
fprintf(fid,'     </coordinates>\n');
fprintf(fid,'   </Point>\n');
fprintf(fid,'  </Placemark>\n');
fprintf(fid,'</Document>\n');
fprintf(fid,'</kml>\n');
 
fclose(fid);
 
return;
