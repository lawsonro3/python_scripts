# This python module contains function to perform common functions

def pv(*vars):
  '''
  This script prints out a vairavle and its value
  '''
  import sys
  try:
    1/0
  except ZeroDivisionError:
    callers_globals = sys.exc_info()[2].tb_frame.f_back.f_globals
    callers_locals = sys.exc_info()[2].tb_frame.f_back.f_locals
  varids=[]  
  i = 0
  for eachvar in vars:
    varids.append(id(eachvar))
  for calling_var_name, calling_var_value in callers_locals.items() +callers_globals.items():
    if id(calling_var_value) in varids:
      if i<1:     
        print calling_var_name, "\t= ", calling_var_value
        i += 1
