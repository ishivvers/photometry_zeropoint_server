'''
Command-line tool providing an interface to the 
 get_SEDs.py library.

'''
from get_SEDs import *
from sys import argv

if __name__ == '__main__':
    filepath = argv[1]
    zp = zeropoint( filepath )
    
    print zp