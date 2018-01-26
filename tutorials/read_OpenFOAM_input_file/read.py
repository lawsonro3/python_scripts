#! /Users/mlawson/env/bin/python3

from python_scripts.openfoam.read import read_input as read

data = read('setUp')
print(data.data)
print(data.file)
