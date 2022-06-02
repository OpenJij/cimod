import os
import sys
import pprint

def test_print_cwd():
  path = os.getcwd()
  pprint.pprint(path)

def test_print_path():
  pprint.pprint(sys.path)
