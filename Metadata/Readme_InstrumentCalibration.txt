'InstrumentCalibration.xlsx' -- gives information about calibration by unique InstrumentID

1. "InstrumentType": Gives type of instrument, which determines method of processing data
2. "InstrumentID”: Unique identifier for instrument. This is preferably the serial number for the instrument, if available, but if not is some other traceable unique identifier.
3. "VarNameGeneral": General variable name for this calibration
4. "CalibrationFactor”: Equals 1 when no calibration is required.  If value differs from one, then multiply by this value when performing calibration.
5. “CalibrationIntercept”: After multiplying by calibration factor, adjust all values by this amount to reach final calibration.