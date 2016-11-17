import sys
import os
import re
import subprocess
import argparse

import settings
import quiver

from kie import KIE_Calculation

import pandas as pd
import glob

def autoquiver(filepath, config_path, gs_p, ts_p, gs_ts_match_p, input_extension='.out', style='g09'):
    if type(gs_p) is str:
        gs_str = gs_p
        gs_p = lambda x: (gs_str in x)
    elif hasattr(gs_p, '__call__'):
        pass
    else:
        raise TypeError("gs_p must either be a string or a function that returns a boolean value.")
    if type(ts_p) is str:
        ts_str = ts_p
        ts_p = lambda x: (ts_str in x)
    elif hasattr(ts_p, '__call__'):
        pass
    else:
        raise TypeError("ts_p must either be a string or a function that returns a boolean value.")
    if type(gs_ts_match_p) is str:
        delimiter = gs_ts_match_p
        gs_ts_match_p = lambda x,y: (x.split(delimiter)[1:] == y.split(delimiter)[1:])
    elif hasattr(ts_p, '__call__'):
        pass
    else:
        raise TypeError("gs_ts_p must either be a string or a function that returns a boolean value.")

    os.chdir(filepath)
    for config in glob.glob("*.config"):
        if os.path.samefile(config_path, config):
            eie_flag = -1
            title = ",,"
            table = ""
            for gs in glob.glob("*"+input_extension):
                if gs_p(gs):
                    for ts in glob.glob("*"+input_extension):
                        if ts_p(ts) and gs_ts_match_p(gs,ts):
                            print "Working on config file: {0}".format(config)
                            print "Using ground state {0} and transition state {1}".format(gs, ts)
                            kie = KIE_Calculation(config, gs, ts, style=style)
                            title_row, row, eie_p = kie.get_row()
                            if eie_flag == -1:
                                eie_flag = eie_p
                            else:
                                if eie_flag != eie_p:
                                    raise ValueError("some calculations represented EIEs and others represented KIEs.")
                            if title is ",,":
                                title = title + title_row
                                table = title + "\n" + table
                            else:
                                if ",," + title_row != title:
                                    raise ValueError("the alignment of the table columns is incorrect.")

                            table += gs + "," + ts + "," + row + "\n"

            with open(os.path.splitext(config)[0]+"-kies.csv", 'w') as f:
                print "----\nAutoQuiver has completed execution.\nOutput is being written to {0}".format(os.path.splitext(config)[0]+"-kies.csv")
                f.write(table)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A program to automatically run PyQuiver on a config file and all ground state and transition states matching certain constraints.")
    parser.add_argument('-v', '--verbose', dest="debug", help='when the verbose flag is set debug information is printed', action='count')
    parser.add_argument('-s', '--style', dest="style", default='g09', help='style of input files')
    parser.add_argument('-e', '--extension', dest="ext", default='.out', help='extension of input files')
    parser.add_argument('config', help='configuration file path')
    parser.add_argument('target', help='target directory file path (where the ground and transition state files live')
    parser.add_argument('gs_p', help='substring in ground state files')
    parser.add_argument('ts_p', help='substring in transition state files')
    parser.add_argument('delimiter', help='delimiter used to match ground and transition state files (all fields separated by the delimiter after the first must match)')

    args = parser.parse_args()
    settings.DEBUG = 0
    if args.debug:
        settings.DEBUG += args.debug
    print settings.DEBUG
        
    autoquiver(args.target, args.config, args.gs_p, args.ts_p, args.delimiter, style=args.style, input_extension=args.ext)
