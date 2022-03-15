#Raw Class has been adapted from the RawFile class in fisher_py. Primarily so that .raw files can be closed after opening them

import numpy as np

from fisher_py.raw_file_reader.raw_file_access import RawFileAccess
from fisher_py.data.filter_enums import MsOrderType, MassAnalyzerType
from fisher_py.data.business import TraceType, ChromatogramTraceSettings, Range, MassOptions
from fisher_py.data import ToleranceUnits, FtAverageOptions, Device

class Raw():
    
    @property
    def first_scan(self) -> int:
        """
        First scan number
        """
        return self._raw_file_access.run_header.first_spectrum
    
    @property
    def last_scan(self) -> int:
        """
        Last scan number
        """
        return self._raw_file_access.run_header.last_spectrum
    
    @property
    def total_time_min(self) -> float:
        """
        Total time of experiment in minutes
        """
        return self._raw_file_access.run_header.end_time
    
    @property
    def number_of_scans(self) -> int:
        """
        Number of scans / spectra
        """
        return self._raw_file_access.run_header_ex.spectra_count
    
    def __init__(self, path):
        self._path = path
        self._raw_file_access = RawFileAccess(path)
        self._raw_file_access.select_instrument(Device.MS, 1)
        
        scan_numbers, rt = self._get_scan_numbers_and_retention_times_()
        self._scan_numbers = scan_numbers
        self._retention_times = rt
        self._spectrum_cache = dict()
        self._result_string_cache = dict()
        
        # fetch retention times and scan numbers for MS1 only
        scan_numbers, rt = self._get_ms_scan_numbers_and_retention_times_(MsOrderType.Ms)
        self._ms1_scan_numbers = scan_numbers
        self._ms1_retention_times = rt
        
        # fetch filters for MS2
        scan_numbers, filter_masses = self._get_ms2_scan_numbers_and_masses_()
        self._ms2_filter_scan_numbers = scan_numbers
        self._ms2_filter_masses = filter_masses
        
    def get_chromatogram(self, mz: float, tolerance: float, trace_type: TraceType=TraceType.MassRange, tolerance_units: ToleranceUnits=ToleranceUnits.ppm, ms_filter: str='ms'):
        """
        Gets chromatogram
        :param mz: Mass/Charge value for mass range chromatogram
        :param tolerace: Tolerance for mass range chromatogram
        :param tolerance_units: Units of the mass tolerance (ppm by default)
        :param ms_filter: Type of MS data (ms or ms2)
        :param trace_type: Type of chromatogram (BasePeek, TIC (total ion current), MassRange (XIC))
        :return: array containing retention times and array containing intensity values
        """

        trace_settings = ChromatogramTraceSettings(trace_type)
        trace_settings.filter = ms_filter
        tolerance_arg = None

        if trace_type == TraceType.MassRange:
            trace_settings.mass_ranges = [Range(mz, mz)]
            tolerance_arg = MassOptions(tolerance, tolerance_units)

        chromatogram_raw = self._raw_file_access.get_chromatogram_data([trace_settings], -1, -1, tolerance_arg)
        return np.array(chromatogram_raw.positions_array[0]), np.array(chromatogram_raw.intensities_array[0])
    
    def _get_scan_(self, scan_number: int):
        mass_analyzer = self._raw_file_access.get_scan_event_for_scan_number(scan_number).mass_analyzer

        if mass_analyzer == MassAnalyzerType.MassAnalyzerFTMS:
            spectrum = self._raw_file_access.get_centroid_stream(scan_number, False)
            positions = np.array(spectrum.masses)
            intensities = np.array(spectrum.intensities)
            charges = np.array(spectrum.charges)
        else:
            spectrum = self._raw_file_access.get_segmented_scan_from_scan_number(scan_number, None)
            positions = np.array(spectrum.masses)
            intensities = np.array(spectrum.intensities)
            charges = np.zeros(positions.shape)
        return positions, intensities, charges
    
    def _get_ms_scan_numbers_and_retention_times_(self, ms_order: MsOrderType):
        first_scan = self.first_scan
        scan_numbers = list()
        retention_times = list()
        scan_events = self._raw_file_access.get_scan_events(self.first_scan, self.last_scan)
        
        for i, scan_event in enumerate(scan_events):
            if scan_event.ms_order != ms_order:
                continue
            scan_number = i + first_scan
            scan_numbers.append(scan_number)
            retention_times.append(self._raw_file_access.retention_time_from_scan_number(scan_number))

        return np.array(scan_numbers), np.array(retention_times)
    
    def _get_ms2_scan_numbers_and_masses_(self):
        first_scan = self.first_scan
        scan_numbers = list()
        filter_mass_values = list()
        scan_events = self._raw_file_access.get_scan_events(self.first_scan, self.last_scan)

        for i, scan_event in enumerate(scan_events):
            if scan_event.ms_order != MsOrderType.Ms2:
                continue
            scan_number = i + first_scan
            precursor_mass = self._get_scan_filter_precursor_mass_(scan_number)
            if precursor_mass is None:
                continue

            scan_numbers.append(scan_number)
            filter_mass_values.append(precursor_mass)

        return scan_numbers, filter_mass_values
    
    def _get_scan_numbers_and_retention_times_(self):
        first_scan = self.first_scan
        scan_numbers = list()
        retention_times = list()

        for i in range(self.number_of_scans):
            scan_number = i + first_scan
            rt = self._raw_file_access.retention_time_from_scan_number(scan_number)
            scan_numbers.append(scan_number)
            retention_times.append(rt)

        return scan_numbers, retention_times
    
    def get_ms1_scan_number_from_retention_time(self, rt: float):
        """
        Get the closest scan number in MS1 spectra from a given retention time (in minutes). If the retention time is smaller than zero
        this will always return the first scan number and if it is larger than the latest entry it will return the biggest
        scan number.
        :param rt: Retention time (in minutes)
        :returns: Scan number
        """
        if rt < 0 or rt > self.total_time_min:
            raise ValueError(f'The retiontion time {rt} is out of bounds. Valid range 0 - {self.total_time_min}.')

        idx = np.argmin(np.abs(self._ms1_retention_times - rt))
        found_rt = self._ms1_retention_times[idx]
        found_scan_nr = self._ms1_scan_numbers[idx]
        return int(found_scan_nr), found_rt
    
    def get_scan_from_scan_number(self, scan_number: int):
        """
        Get scan data from a scan number. The data returned is structured in a tuple as follows:
            (masses, intensities, ion_charges, scan_event_descriptions)
        :param scan_number: The number of the scan
        :returns: Tuple organized as (masses, intensities, ion_charges, scan_event_descriptions)
        """
        if scan_number < self.first_scan or scan_number > self.last_scan:
            raise ValueError(f'The scan number {scan_number} is out of bounds. Valid range {self.first_scan} - {self.last_scan}.')

        if scan_number in self._spectrum_cache:
            return self._spectrum_cache[scan_number][0], self._spectrum_cache[scan_number][1], self._spectrum_cache[scan_number][2], self._result_string_cache[scan_number]

        positions, intensities, charges = self._get_scan_(scan_number)
        scan_event_str = self._raw_file_access.get_scan_event_string_for_scan_number(scan_number)
        return positions, intensities, charges, scan_event_str
    
    def get_scan_ms1(self, rt: float):
        """
        Gets MS1 (MS) spectrum for a given retention time value in minutes
        :param rt: retention time value in minutes
        :returns: Three arrays containing Mass/Charge values, intensity values and charge values as well as the actual retention time
        """
        scan_number, found_rt = self.get_ms1_scan_number_from_retention_time(rt)
        rts, intensities, charges, _ = self.get_scan_from_scan_number(scan_number)
        return rts, intensities, charges, found_rt
    
    def get_scan_ms2(self, rt: float, precursor_mz: float=None):
        """
        Gets MS2 spectrum for a given retention time value in minutes
        :param rt: retention time value in minutes
        :param recursor_mz: Optional precursor mass value to filter spectra
        :returns: Three arrays containing Mass/Charge values, intensity values and charge values as well as the actual retention time
        """
        scan_number, found_rt = self.get_ms2_scan_number_from_retention_time(rt, precursor_mz)
        rts, intensities, charges, _ = self.get_scan_from_scan_number(scan_number)
        return rts, intensities, charges, found_rt
    
    def _get_scan_filter_precursor_mass_(self, scan_number: int):
        scan_event = self._raw_file_access.get_scan_event_for_scan_number(scan_number)
        if scan_event.ms_order != MsOrderType.Ms2:
            return None
        return scan_event.get_reaction(0).precursor_mass
    
    def get_ms2_scan_number_from_retention_time(self, rt: float, precursor_mz: float=None, tolerance_ppm = 10e-3):
        """
        Get the closest scan number in MS2 spectra from a given retention time (in minutes). If the retention time is smaller than zero
        this will always return the first scan number and if it is larger than the latest entry it will return the biggest
        scan number.
        NOTE: This method does not yet support all mass tolerance units
        :param rt: Retention time (in minutes)
        :param precursor_mz: Precursor mass
        :param tolerance: Mass tolerance to precursor (in ppm)
        :returns: Scan number
        """
        if rt < 0 or rt > self.total_time_min:
            raise ValueError(f'The retiontion time {rt} is out of bounds. Valid range 0 - {self.total_time_min}.')

        if precursor_mz is None:
            idx = np.argmin(np.abs(self._ms2_retention_times - rt))
            found_rt = self._ms2_retention_times[idx]
            found_scan_nr = int(self._ms2_scan_numbers[idx])
            return found_scan_nr, found_rt
        else:
            ms2_idxs = np.abs(self._ms2_filter_masses - precursor_mz) <= tolerance_ppm
            ms2_scans = self._ms2_filter_scan_numbers[ms2_idxs]
            ms2_rts = self._retention_times[ms2_scans.astype(int) + self.first_scan]
            idx = np.argmin(np.abs(ms2_rts - rt))
            found_rt = ms2_rts[idx]
            found_scan_nr = int(ms2_scans[idx])
            return found_scan_nr, found_rt
			
    def close(self):
        self._raw_file_access.dispose()
    