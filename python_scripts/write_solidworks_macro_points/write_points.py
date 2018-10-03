import os
import numpy as np

def wirte_macro(input_file,output_file):

    with open(input_file,'r') as fid:
            raw = fid.readlines()

    with open(output_file,"a") as write_output:
        write_output.write('Part.InsertCurveFileBegin\n')
        for i, line in enumerate(raw):
            write_output.write(''.join(['boolstatus = Part.InsertCurveFilePoint(',line.split()[0],',',line.split()[1],',',line.split()[2],')\n']))
        write_output.write('boolstatus = Part.InsertCurveFileEnd()\n')
        write_output.write(''.join(['boolstatus = Part.SelectedFeatureProperties(0, 0, 0, 0, 0, 0, 0, 1, 0, "',input_file[-10:-7],'")\n']))
        write_output.write('Part.ClearSelection2 True\n')

if __name__ == '__main__':
    loc = '/Users/mlawson/GoogleDrive/Work/NREL/Projects/HFM-ECP/nrel_5mw/airfoil-files-prssure-suction/'
    output_file = loc + '2_10_sw_macro.txt'
    header = open(output_file,'w+')
    header.write('\n')
    header.write('Dim swApp As Object\n')
    header.write('\n')
    header.write('Dim Part As Object\n')
    header.write('Dim boolstatus As Boolean\n')
    header.write('Dim longstatus As Long, longwarnings As Long\n')
    header.write('\n')
    header.write('Sub main()\n')
    header.write('\n')
    header.write('Set swApp = Application.SldWorks\n')
    header.write('\n')
    header.write('Set Part = swApp.ActiveDoc\n')
    header.write('Dim myModelView As Object\n')
    header.write('Set myModelView = Part.ActiveView\n')
    header.write('myModelView.FrameState = swWindowState_e.swWindowMaximized\n')
    header.close()

    aero_file = '/Users/mlawson/Documents/GitHub/python_scripts/python_scripts/nrel_5mw/aero_properties.yaml'
    tmp = open('aero_file','r')
    import yaml
    aero=yaml.load(tmp)
    tmp.close()



    for i in range(2,11):
        loc = '/Users/mlawson/GoogleDrive/Work/NREL/Projects/HFM-ECP/nrel_5mw/airfoil-files-prssure-suction/'
        input_file = ''.join([loc,'r',str(i),'_P.sldcrv'])
        wirte_macro(input_file,output_file,z)

    footer = open(output_file,'a')
    footer.write('Sub End')
