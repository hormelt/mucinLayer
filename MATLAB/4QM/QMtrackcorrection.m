function QMTracks = QMtrackcorrection(Tracks,bpData,refCenters,CalibParams, ...
                                      FeatSize,DeltaFit)

% Provides corrections for a pretrack based on the 4QM method.
%
% INPUT:
%   Tracks: The tracks obtained by classical tracking.
%   bpData: The bandpassed image data.
%   refCenters: The refCenters of all the particles.
%   CalibParams: The CalbrationParameters obtained by calibration.
%   NParticles: Number of Particles.
%   FeatSize: The Size of the particles.
%   DeltaFit: A value to adjust the size of the subframe around the
%   particle.
%
% OUTPUT:
%   QMTracks: The Tracks corrected using the 4QM method.


QMTracks = zeros(size(Tracks,1),5);
QMTrackRow = 1; % for indexing
for i = 1:size(CalibParams,1)
    ParticleID = CalibParams(i,7);
    
    % Construct subdata
    TrackFrames = Tracks(Tracks(:,6)==ParticleID,5);
    NTrackFrames = numel(TrackFrames);
    xCoarse = refCenters(ParticleID,1);
    yCoarse = refCenters(ParticleID,2);      
    Cols = SetAxisSubdata(xCoarse,FeatSize,DeltaFit);
    Rows = SetAxisSubdata(yCoarse,FeatSize,DeltaFit);      
    subData = bpData(Rows,Cols,TrackFrames);
    
    % Calculate the true centers of the particles.
    pCoef = CalibParams(i,:);
    [A,B,C,D] = FQM(subData);
    TrackCorrection = [pCoef(1)*(-A-C+B+D)./(A+B+C+D)+pCoef(2) ...
                       pCoef(4)*(-A-B+C+D)./(A+B+C+D)+pCoef(5) ...
                       (1:size(subData,3))'];
    preTrack = [xCoarse*ones(NTrackFrames,1), ...
                yCoarse*ones(NTrackFrames,1), ...
                zeros(NTrackFrames,1)];
    correctedTrack = preTrack + TrackCorrection;
    QMTracks(QMTrackRow:QMTrackRow+NTrackFrames-1,:) = [correctedTrack ParticleID*ones(numel(TrackFrames),1),...
        CalibParams(i,8)*ones(numel(TrackFrames),1)];
    QMTrackRow = QMTrackRow+NTrackFrames;
end



