import numpy as np
from ctypes import CDLL, POINTER, cast, c_uint, c_double, c_int, c_char_p
import pathlib
spclib = CDLL(str(pathlib.Path(__file__).parent.resolve()/"spclib.so"))
doublePtr = POINTER(c_double)
doublePtrPtr = POINTER(doublePtr)
uintPtr = POINTER(c_uint)
spclib.spc.argtypes = [doublePtrPtr,c_char_p, uintPtr, uintPtr]
spclib.spc.restype = c_int



class SPC():
    def __init__( self,
        mintemp=0, maxtemp=0.251,  tempstep=0.01,
        swcycles=100, nearest_neighbours=11,
        mstree=True, directedgrowth=True,
        ncl_reported=12,
        randomseed=0
    ):
        assert mintemp<=maxtemp, "error: mintemp > maxtemp" 
        self.mintemp = mintemp
        self.maxtemp = maxtemp
        self.tempstep = tempstep
        self.swcycles = swcycles
        self.nearest_neighbours = nearest_neighbours
        self.ncl_reported = ncl_reported
        self.mstree_char = _bool2char(mstree)
        self.directedgrowth_char = _bool2char(directedgrowth)
        self.randomseed = randomseed

        self.temp_vector = np.arange(self.mintemp,self.maxtemp,self.tempstep)

        #all additional ouputs disabled
        self.output_name = 'spc_data'
        self.writelabels_char = _bool2char(False)
        self.writecorfile_char = _bool2char(False)
        self.savesuscept_char = _bool2char(False)
        self.writefpsum_char = _bool2char(False)
        self.writesizes_char = _bool2char(False)

    @property
    def temp_vectr(self):
        return self.temp_vector
    
    def config_additional_outputs(self, output_name, writelabels,  writecorfile, 
        savesuscept, writefpsum, writesizes):
        self.output_name = output_name,
        self.writelabels_char = _bool2char(writelabels)
        self.writecorfile_char = _bool2char(writecorfile)
        self.savesuscept_char = _bool2char(savesuscept)
        self.writefpsum_char = _bool2char(writefpsum)
        self.writesizes_char = _bool2char(writesizes)
        
    def run(self, data, return_sizes=False):
        dim = data.shape[1]
        data_points = data.shape[0]
        
        ntemp = self.temp_vector.shape[0]

        classes = np.ascontiguousarray(np.zeros((ntemp,data_points),dtype=c_uint, order="C"))

        new_param = (f"NumberOfPoints: {data_points}\n"
                f"OutFile: {self.output_name}\n"
                f"Dimensions: {dim}\n"
                f'MinTemp: {self.mintemp}\n'
                f'MaxTemp: {self.maxtemp}\n'
                f'TempStep: {self.tempstep}\n'
                f'SWCycles: {self.swcycles}\n'
                f'KNearestNeighbours: {self.nearest_neighbours}\n'
                f'MSTree{self.mstree_char}\n'
                f'DirectedGrowth{self.directedgrowth_char}\n'
                f'WriteLabels{self.writelabels_char}\n'
                f"ClustersReported: {self.ncl_reported}\n"
                f'WriteCorFile{self.writecorfile_char}\n'
                f"SaveSuscept{self.savesuscept_char}\n"
                f"WriteFPSum{self.writefpsum_char}\n"
                f"WriteSizes{self.writesizes_char}\n"
                f"ForceRandomSeed: {self.randomseed}\n"
                "Timing~\n" 
                ) 
        
        ct_arr = np.ctypeslib.as_ctypes(data.astype(c_double, order='C',copy=False)) 
        
        doublePtrArr = doublePtr * data_points
        input_data = cast(doublePtrArr(*(cast(row, doublePtr) for row in ct_arr)), doublePtrPtr)


        if return_sizes:
            sizes = np.ascontiguousarray(np.zeros((ntemp,self.ncl_reported),dtype=c_uint, order="C"))

        res = spclib.spc(input_data, new_param.encode(),
            cast(np.ctypeslib.as_ctypes(classes),uintPtr ),
            cast(np.ctypeslib.as_ctypes(sizes), uintPtr) if return_sizes else cast(None, uintPtr))
        if res != 0:
            raise(f'spc c code returned: {res}')
        
        if return_sizes:
            return classes, sizes
        return classes

def _bool2char(x):
    return '|' if x else '~'