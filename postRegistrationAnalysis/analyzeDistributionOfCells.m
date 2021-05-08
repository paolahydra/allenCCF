% 622	'Bed nuclei of the stria terminalis'
% 1260	'stria terminalis'
% 572	'Striatum'

% 598	'Central amygdalar nucleus lateral part'
% 599	'Central amygdalar nucleus medial part'
% 597	'Central amygdalar nucleus capsular part'


% 801   'Subthalamic nucleus'
% 798   'Parasubthalamic nucleus'
% 1222  'nigrostriatal tract'


% 795	'Lateral hypothalamic area'
% 716	'Hypothalamus'
% 793	'Posterior hypothalamic nucleus'
% 721	'Paraventricular hypothalamic nucleus'

% 803	'Zona incerta'
% 805	'Fields of Forel'
% 646	'Ventral medial nucleus of the thalamus'


% 823	'Substantia nigra reticular part'
% 867	'Substantia nigra compact part'

% 826	'Midbrain reticular nucleus retrorubral area'
% 827	'Midbrain reticular nucleus'
% 807	'Midbrain'

% 853	'Cuneiform nucleus'
% 868	'Pedunculopontine nucleus'


%% load data to analyze
folders2Include = uipickfiles('FilterSpec', '/Users/galileo/dati/registered_brains_completed');
plotFWire = 1;

S = loadTabDataFromMultipleBrains(folders2Include, plotFWire); %update to exclude point that are in 'root' (avIndex = 1)

%% find unique targets and sort them
i = 1;


all_Names = S(i).T_roi.name;
all_Files = S(i).T_roi.roiFIle;
all_acIndex = S(i).T_roi.avIndex;
all_isPositiveML_location = S(i).T_roi.ML_location;
all_isPositiveML_location = all_isPositiveML_location>0; %this should be on the right

TargetAreas_R = unique(all_acIndex(all_isPositiveML_location));
for n = 1:length(TargetAreas_R)
    TargetNames_R{n,1} = all_Names{find(all_acIndex==TargetAreas_R(n), 1)};
    % count how many registered slices contributed to any given area count
    TargetSliceCount_R{n,1} = length(unique(all_Files(all_acIndex==TargetAreas_R(n) & all_isPositiveML_location)));
end

TargetAreas_L = unique(all_acIndex(~all_isPositiveML_location));
for n = 1:length(TargetAreas_L)
    TargetNames_L{n,1} = all_Names{find(all_acIndex==TargetAreas_L(n), 1)};
    % count how many registered slices contributed to any given area count
    TargetSliceCount_L{n,1} = length(unique(all_Files(all_acIndex==TargetAreas_L(n) & ~all_isPositiveML_location)));
end

%%

tbl = tabulate(all_acIndex(all_isPositiveML_location));
tbl(tbl(:,2)==0,:) = [];
[~, i_sorting] = sort(tbl(:,2), 'descend');
tbl = tbl(i_sorting,:);
TargetNames_R = TargetNames_R(i_sorting, :);
TargetSliceCount_R = TargetSliceCount_R(i_sorting, :);
S(i).cellCounts_right = table(tbl(:,1), TargetNames, tbl(:,2), TargetSliceCount, 'VariableNames', {'avIndex','name','count','nSlicesDenominator'});

tbl = tabulate(all_acIndex(~all_isPositiveML_location));
tbl(tbl(:,2)==0,:) = [];
[~, i_sorting] = sort(tbl(:,2), 'descend');
tbl = tbl(i_sorting,:);
TargetNames = TargetNames(i_sorting, :);
S(i).cellCounts_left = table(tbl(:,1), TargetNames, tbl(:,2), TargetSliceCount, 'VariableNames', {'avIndex','name','count','nSlicesDenominator'});


% for a given region of interest that I want to compare, I need to reload
% each registered slice that contained that area, mask the area, and sum
% the extent of pixels across all slices (as a proxy of the volume
% considered).
%
% for now, the number of slices will do to have a very coarse sense (but
% it's actuall not even good because I would need the number of slices that
% intersect that region of interest, not the number of slices that have at
% least one cell
