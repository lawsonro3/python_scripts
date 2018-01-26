import os

import numpy as np

class OpenFOAMInput(object):
    '''OpenFOAM input data object
    '''
    def __init__(self,input_file):
        self.data = {}
        self.file = os.path.abspath(input_file)
        self._read()


    def _read(self):
        '''
        Note that this function only reads scalar and boolean/string inputs from
        an OpenFOAM intput input file
        '''
        
        with open(self.file,'r') as fid:
            raw = fid.readlines()

        count = 0
        bloc_comment_test = False
        for i, line in enumerate(raw):

            if raw[i][0:2]=='/*':
                bloc_comment_test = True

            if bloc_comment_test is False:

                if raw[i].strip()[0:2]=='//' or raw[i].strip()[0:1]=='#': # Check if the string is a comment and skip line
                    pass

                elif len(raw[i].strip())==0: # Check if the string is empty and skip line
                    pass

                else:
                    tmp = raw[i].strip().rstrip().split()
                    try:
                        self.data[tmp[0]] = np.float(tmp[1][:-1])
                    except:
                        self.data[tmp[0]] = tmp[1][:-1]



            if raw[i][0:2]=='\*':
                bloc_comment_test = False


def read_input(input_file):
    print('\nReading ' + input_file)

    data = OpenFOAMInput(input_file)

    return data
