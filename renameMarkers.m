function [newMarkers] = renameMarkers(markers)
    
        %Left side markers
        try
            newMarkers.LPSI = markers.L_PSI;
        catch 
            newMarkers.LPSI = markers.L_IPS;
        end
        newMarkers.LASI = markers.L_ASI;
        newMarkers.LFM1 = markers.L_FM1;
        newMarkers.LFM2 = markers.L_FM2;
        newMarkers.LFM5 = markers.L_FM5;
        newMarkers.LFCC = markers.L_FCC;
        newMarkers.LTAM = markers.L_TAM;
        newMarkers.LFAL = markers.L_FAL;
        newMarkers.LRSP = markers.L_RSP;
        newMarkers.LUSP = markers.L_USP;

        % Right side markers
        try
            newMarkers.RPSI = markers.R_PSI;
        catch
            newMarkers.RPSI = markers.R_IPS;
        end
        newMarkers.RASI = markers.R_ASI;
        newMarkers.RFM1 = markers.R_FM1;
        newMarkers.RFM2 = markers.R_FM2;
        newMarkers.RFM5 = markers.R_FM5;
        newMarkers.RFCC = markers.R_FCC;
        newMarkers.RTAM = markers.R_TAM;
        newMarkers.RFAL = markers.R_FAL;
        newMarkers.RRSP = markers.R_RSP;
        newMarkers.RUSP = markers.R_USP;

        % Trunk
        newMarkers.CLAV = markers.CLAV;
        newMarkers.TV10 = markers.TV10;

        % Shoulder
        newMarkers.RSHO = markers.R_SHO;
        newMarkers.LSHO = markers.L_SHO;
        
end