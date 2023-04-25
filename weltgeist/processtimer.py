"""
Make an object that times processes
Sam Geen, April 2023
"""

from time import process_time

class ProcessTimer():
    def __init__(self):
        self._totals = {}
        self._begins = {}
        self._namestack = ""

    def Begin(self,name):
        """
        Begin the time checking for named process

        Parameters
        ----------
        name: string
            name of process to time
        """
        self._namestack = self._namestack+"."+name
        self._begins[self._namestack] = process_time()


    def End(self,name):
        """
        End the time checking for named process

        Parameters
        ----------
        name: string
            name of process to time
        """
        end = process_time()
        namestack = self._namestack
        diff = end - self._begins[namestack]
        if not namestack in self._totals:
            self._totals[namestack] = 0.0
        self._totals[namestack] += diff
        self._namestack = namestack[:-len("."+name)]

    def OutputLog(self,filename):
        """
        Write the totals in the timer to a log file

        Parameters
        ----------
        filename: string
            name of text file to output to
        """
        f = open(filename,"w")
        lcol = 30 # Length of column
        f.write("Process name".ljust(lcol)+" : Time in seconds\n")
        # TODO: sort somehow so the stack is well ordered
        for key, value in self._totals.items():
            f.write(key.ljust(lcol)+" : "+str(value)+"\n")
        f.close()

    def Save(self, h5file):
        """
        Save internal state 
        """
        # TODO: this
        # Suggest refactoring how integrator saves to include a wrapper around h5py
        raise NotImplementedError

    def Load(self, h5file):
        """
        Load internal state 
        """
        # TODO: this
        # Suggest refactoring how integrator saves to include a wrapper around h5py
        raise NotImplementedError
