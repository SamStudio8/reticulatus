from datetime import datetime
import sys
import os

FLAG_PATH = sys.argv[1]
for line in open(os.path.join(FLAG_PATH, "q.txt")):
    fields = line.strip().split()
    sample = fields[0]
    depth = fields[1]
    status = fields[2]
    if status == '-':
        continue
    jobid = fields[3]
    rule = fields[4]

    start_p = os.path.join(FLAG_PATH, '%s.start' % jobid)
    finish_p = os.path.join(FLAG_PATH, '%s.finish' % jobid)

    dts = [None, None]
    times = ['-', '-']
    if os.path.exists(start_p):
        dts[0] = datetime.fromtimestamp(os.path.getmtime(start_p))
        times[0] = dts[0].strftime('%Y-%m-%d_%H:%M:%S')

    if os.path.exists(finish_p):
        dts[1] = datetime.fromtimestamp(os.path.getmtime(finish_p))
        times[1] = dts[1].strftime('%Y-%m-%d_%H:%M:%S')

    if dts[0] and dts[1] and dts[0] > dts[1]:
        times[1] = '-'

    fields.extend(times)
    print("\t".join(fields))
