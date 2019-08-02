import os
import datetime
import numpy as np

import InstrumentDriver


class Driver(InstrumentDriver.InstrumentWorker):
    """This class reads the last entry in the BlueFors log files."""

    def _readLastLine(self, filename_prefix):
        """Returns the last line in the log file.

        Parameters
        ----------
        filename_prefix : str
            Filename prefix corresponding to the physical quantity.

        Returns
        -------
        line : str
            Last line in the log file.

        Raises
        ------
        FileNotFoundError
            No log file found.
        """
        # Get the list of files in the folder and return the one with
        # the most recent name.
        subfolders = sorted([f for f in os.listdir(self._path)
                    if not os.path.isfile(os.path.join(self._path, f))])

        # Check the recent folders for the log file.
        not_found = True
        for i in range(-1, -len(subfolders) - 1, -1):
            filename = ''.join([filename_prefix, subfolders[i], '.log'])
            full_filename = os.path.join(self._path, subfolders[i], filename)
            if os.path.isfile(full_filename):
                not_found = False
                break
        if not_found:
            return None

        # Read the last line in the log file.
        offset = 1024
        with open(full_filename, 'rb') as f:
            f.seek(0, os.SEEK_END)
            # Get the size of the file.
            sz = f.tell()
            while True:
                if offset > sz:
                    offset = sz
                f.seek(-offset, os.SEEK_END)
                lines = f.readlines()
                if len(lines) > 1 or offset == sz:
                    line = lines[-1]
                    break
                offset *= 2

        # Compare the record timestamp with the current time.
        str_time = b','.join(line.split(b',', 2)[:2]).decode().strip()
        dt_time = datetime.datetime.strptime(str_time,
                                             '%y-%m-%d,%H:%M:%S')
        now = datetime.datetime.now()
        time_lapsed = (now - dt_time).total_seconds()
        if time_lapsed > 120:
            return None
        return line

    def _parseSimpleRecord(self, line=None):
        """Parses a thermometer/flowmeter record.

        Parameters
        ----------
        line : None or str
            Single thermometer log record. If None (default), the method
            will return np.NaN.

        Returns
        -------
        float
            Extracted value.
        """
        if line is None:
            return np.NaN
        return float(line.split(b',')[-1])

    def _parsePressureRecord(self, line=None, channel=1, field='status'):
        """Parses a pressure record.

        Parameters
        ----------
        line : str or None
            Single thermometer log record. If None (default), the method
            will return np.NaN.
        channel : {1, 2, 3, 5, 6}
            Maxigauge channel number.
        field: {'status', 'pressure'}
            Value to return, i.e., either the status of the gauge or
            its pressure reading.

        Returns
        -------
        float or bool
            Extracted value of the field.
        """
        if line is None:
            return np.NaN
        fields = line.split(b',')
        if field == 'status':
            return bool(int(fields[6 * channel - 2].decode()))
        elif field == 'pressure':
            # Extract the pressure and convert it from mbar to bar.
            if self._parsePressureRecord(line, channel, 'status'):
                return float(fields[6 * channel - 1].decode()) / 1e3
            else:
                return np.NaN
        else:
            raise ValueError("Parameter 'field' could take only "
                             "'status' or 'pressure' value.")

    def performGetValue(self, quant, options={}):
        """Perform the Get Value for the driver operation."""
        self._path = self.getValue('Path to BlueFors Logs')
        if quant.name == '50K Flange Temperature':
            # Read the last line in the log file.
            line = self._readLastLine('CH1 T ')
            # Extract the temperature.
            return self._parseSimpleRecord(line)
        elif quant.name == '50K Flange Thermometer Resistance':
            # Read the last line in the log file.
            line = self._readLastLine('CH1 R ')
            # Extract the resistance.
            return self._parseSimpleRecord(line)
        elif quant.name == '4K Flange Temperature':
            # Read the last line in the log file.
            line = self._readLastLine('CH2 T ')
            # Extract the temperature.
            return self._parseSimpleRecord(line)
        elif quant.name == '4K Flange Thermometer Resistance':
            # Read the last line in the log file.
            line = self._readLastLine('CH2 R ')
            # Extract the resistance.
            return self._parseSimpleRecord(line)
        elif quant.name == 'Still Temperature':
            # Read the last line in the log file.
            line = self._readLastLine('CH5 T ')
            # Extract the temperature.
            return self._parseSimpleRecord(line)
        elif quant.name == 'Still Thermometer Resistance':
            # Read the last line in the log file.
            line = self._readLastLine('CH5 R ')
            # Extract the resistance.
            return self._parseSimpleRecord(line)
        elif quant.name == 'Mix Temperature':
            # Read the last line in the log file.
            line = self._readLastLine('CH6 T ')
            # Extract the temperature.
            return self._parseSimpleRecord(line)
        elif quant.name == 'Mix Thermometer Resistance':
            # Read the last line in the log file.
            line = self._readLastLine('CH6 R ')
            # Extract the resistance.
            return self._parseSimpleRecord(line)
        elif quant.name in ('P1 Gauge Status', 'P2 Gauge Status',
                            'P3 Gauge Status', 'P4 Gauge Status',
                            'P5 Gauge Status', 'P6 Gauge Status'):
            # Read the last line in the log file.
            line = self._readLastLine('maxigauge ')
            # Extract the gauge status.
            channel = int(quant.name[1])
            return self._parsePressureRecord(line, channel, 'status')
        elif quant.name in ('P1 Pressure', 'P2 Pressure',
                            'P3 Pressure', 'P4 Pressure',
                            'P5 Pressure', 'P6 Pressure'):
            # Read the last line in the log file.
            line = self._readLastLine('maxigauge ')
            # Extract the pressure.
            channel = int(quant.name[1])
            return self._parsePressureRecord(line, channel, 'pressure')
        elif quant.name == 'Flow Rate':
            # Read the last line in the log file.
            line = self._readLastLine('Flowmeter ')
            # Extract the flow rate and convert it to mol/s from
            # mmol/s.
            return self._parseSimpleRecord(line) / 1.e3
        else:
            return quant.getValue()


if __name__ == '__main__':
    pass
