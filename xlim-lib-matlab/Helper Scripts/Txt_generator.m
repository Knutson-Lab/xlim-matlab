binWidth = input("Binwidth: ");
channels = input("Channels: ");
NumberIRF = input("Number of irf curves: ");
NumberData = input("Number of data curves: ");
NumberTau = input("Number of tau values: ");
tau = zeros(NumberTau,1);
Photons = zeros(NumberTau,1);
for i = 1:NumberTau
    tauprompt = append('The value of tau',int2str(i), ': ');
    tau(i,1) = input(tauprompt);
    photonprompt = append('The value of photon',int2str(i), ': ');
    Photons(i,1) = input(photonprompt);
end

Filename = "Simulation_Parameters.txt";
fileID = fopen(Filename,'w');
fprintf(fileID,'Binwidth: %f\n', binWidth);
fprintf(fileID,'Channels: %d\n', channels);
fprintf(fileID,'Number of irf curves: %d\n', NumberIRF);
fprintf(fileID,'Number of data curves: %d\n', NumberData);
fprintf(fileID,'Number of tau values: %d\n', NumberData);
for i = 1:NumberTau
    tauprompt = append('The value of tau',int2str(i), ': %f\n');
    fprintf(fileID,tauprompt, tau(i));
    photonprompt = append('The value of photon',int2str(i), ': %d\n');
    fprintf(fileID,photonprompt, Photons(i));
end
fclose(fileID);
irf = sim_tcspc_irfs(binWidth,zeros(channels,1),[20000, 0.8, 0.25, 0, 0],NumberIRF);
switch NumberTau
    case 1
        datacurves = sim_tcspc_dks(binWidth,irf,zeros(channels,1), ...
            [Photons(1,1), tau(1,1), 0, 0, 0, 0],NumberTau,NumberData);
    case 2
        datacurves = sim_tcspc_dks(binWidth,irf,zeros(channels,1), ...
            [Photons(1,1), tau(1,1), Photons(2,1), tau(2,1),0, 0, ...
            0, 0],NumberTau,NumberData);
    case 3
        datacurves = sim_tcspc_dks(binWidth,irf,zeros(channels,1), ...
            [Photons(1,1), tau(1,1), Photons(2,1), tau(2,1), ...
            Photons(3,1), tau(3,1), 0, 0, 0, 0],NumberTau,NumberData);
    case 4
        datacurves = sim_tcspc_dks(binWidth,irf,zeros(channels,1), ...
            [Photons(1,1), tau(1,1), Photons(2,1), tau(2,1), ...
            Photons(3,1), tau(3,1), Photons(4,1), tau(4,1), ...
            0, 0, 0, 0],NumberTau,NumberData);
    case 5
        datacurves = sim_tcspc_dks(binWidth,irf,zeros(channels,1), ...
            [Photons(1,1), tau(1,1), Photons(2,1), tau(2,1), ...
            Photons(3,1), tau(3,1), Photons(4,1), tau(4,1), ...
            Photons(5,1), tau(5,1),0, 0, 0, 0],NumberTau,NumberData); 
end
SizeV = size(irf);
for i = 1:SizeV(2)
    Filename = append('irf',int2str(i),'.txt');
    fileID = fopen(Filename,'w');
    fprintf(fileID,'%d\n',channels);
    fprintf(fileID,'%d\n',irf(:,i));
    fclose(fileID);
end
SizeV = size(datacurves);
for i = 1:SizeV(2)
    Filename = append('data_curves',int2str(i),'.txt');
    fileID = fopen(Filename,'w');
    fprintf(fileID,'%d\n',channels);
    fprintf(fileID,'%d\n',datacurves(:,i));
    fclose(fileID);
end
