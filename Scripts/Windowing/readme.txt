-- SCRIPTS --

CreateTimeWindows.m (Windowing) => TimeWindows_30min.mat (ProcessedData) - denotes start times, end times, and dates for windows

CreateDataWindows.m (Windowing) => DataWindows_30min.mat (ProcessedData) - breaks wind and flux data into 30 minute windows

CreateWindowAverage.m (Windowing) => WindowAverageWindows.mat (ProcessedData) - creates window average time series at multiple time scales for 30 minute windows

CreateActivityWindows.m (Windowing) => ActivityWindows.mat - determines flux activities (fQ) from window average windows, as well as corresponding TFEM thresholds and hysteresis analyses

CreateLowPass.m (Windowing) => LowPassWindows.mat (ProcessedData) - creates low pass time series at multiple time scales for 30 minute windows

CreateThresholdWindows.m (Windowing) - Create all data for threshold analysis - OUT OF DATE!

CreateFluxLawWindows.m (Windowing) - Create all data for flux law analysis - OUT OF DATE!

-- FUNCTIONS --

IntersectingTimeIntervals.m

CreateTimeBlocks.m

LowPassFilter.m

-- OTHER --

"Old" - past scripts no longer in use