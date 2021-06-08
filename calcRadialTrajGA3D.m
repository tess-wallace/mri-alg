function [trajRAD_mm, trajRAD, polarAngle, azimuthalAngle ] = calcRadialTrajGA3D(params)


% Compute golden means ratios in 3D
Mfib3d = [0, 1, 0; 0, 0, 1; 1, 0, 1];

[V3d, D3d] = eig(Mfib3d);

v = V3d(:,1)/V3d(3,1);

m_phi1 = round(v(1), 15);
m_phi2 = round(v(2), 15);

baseResolution = params.baseResolution;
nCol = params.nCol;
nSpokes = params.nSpokes;
fovx_mm = params.fovx_mm;
osFactor = nCol/baseResolution;

maxkspace=1/(fovx_mm/baseResolution)/2;

inc = 2*maxkspace/(baseResolution*osFactor); 

rho =-maxkspace(1):inc(1):maxkspace(1)-inc(1); rho=rho.';

% kspace trajectory
trajx = zeros(nCol,nSpokes);
trajy = zeros(nCol,nSpokes);
trajz = zeros(nCol,nSpokes);

polarAngle = zeros(1,nSpokes);

azimuthalAngle = zeros(1,nSpokes);

for ii=1:nSpokes
    
    if ii == 1
        m1(ii) = 0;
        m2(ii) = 0;
    else
        m1(ii) = mod( (m1(ii-1) + m_phi1), 1);
        m2(ii) = mod( (m2(ii-1) + m_phi2), 1);
    end
    
    polarAngle(ii) = pi/2 + m1(ii)*pi/2;
    
    azimuthalAngle(ii) = 2*pi*m2(ii);

    xA = cos(azimuthalAngle(ii))*sin(polarAngle(ii));
    
    yA = sin(azimuthalAngle(ii))*sin(polarAngle(ii));
    
    zA = cos(polarAngle(ii));
    
    trajx(:,ii)= rho*xA;

    trajy(:,ii)= rho*yA;

    trajz(:,ii)= rho*zA;
        	
end

trajRAD(:,1) = reshape(trajx,[],1);
trajRAD(:,2) = reshape(trajy,[],1);
trajRAD(:,3) = reshape(trajz,[],1);

trajRAD_mm = trajRAD; % original

trajRAD=trajRAD/(2*max(max(abs(trajRAD)))); % scaled between -0.5 and 0.5