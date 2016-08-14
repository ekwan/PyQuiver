import sys
import os
import re
import subprocess

from kie import KIE_Calculation

import pandas as pd
import glob

def auto_quiver(filepath, gs_p, ts_p, gs_ts_match_p, input_extension='.snip', style='g09'):
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
        eie_flag = -1
        title = ",,"
        table = ""
        for gs in glob.glob("*"+input_extension):
            if gs_p(gs):
                for ts in glob.glob("*"+input_extension):
                    if ts_p(ts) and gs_ts_match_p(gs,ts):
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
            f.write(table)

if __name__ == "__main__":
    auto_quiver(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
