import numpy as np

def read_log(log_file=None):
    '''This function reads Nalu log files
       Currently, the function only reads timing info output by nalu-wind
       It would be good to add more functionality to this function

    '''

    if log_file is None:
        raise Exception('Please enter a log file name')

    with open(log_file,'r') as fid:
        raw = fid.readlines()

    count = 0
    for i, line in enumerate(raw):

        # Read timing information from the log file
        if np.size(line.split()) == 0:
            pass
        elif line.split()[0]  == 'WallClockTime:':
            tmp = line.split()
            if count == 0:
                time_headers = [tmp[0],tmp[2],tmp[4],tmp[6],tmp[8]]
                times = np.array([[tmp[1],tmp[3],tmp[5],tmp[7],tmp[9]]])
            else:
                times = np.append(times,[[tmp[1],tmp[3],tmp[5],tmp[7],tmp[9]]],axis=0)

            count += 1

    times = times.astype(np.float)
    return time_headers,times
