import os
import sys

LOG_DIR = sys.argv[1]

search = {
    "loaded sequences": {'f': lambda x, struct: x.strip().split()[-2]},
    "loaded overlaps": {'f': lambda x, struct: x.strip().split()[-2]},
    "aligning overlaps": {'f': lambda x, struct: x.strip().split()[-2]},
    "generating consensus": {'f': lambda x, struct: x.strip().split()[-2]},
    "polished windows on GPU": {'f': lambda x, struct: x.strip().split()[-2]},
    "polished remaining windows on CPU": {'f': lambda x, struct: x.strip().split()[-2]},
    "generated consensus": {'f': lambda x, struct: x.strip().split()[-2]},
    "total =": {"name": "total", 'f': lambda x, struct: x.strip().split()[-2]},
    #"With cudnn": lambda x: x.strip().split()[-1],
    #"Setting tensorflow threads to": lambda x: x.strip().split()[-1],
    #"Predict] Processing region": lambda x: x.strip().split()[0].replace('[',''),
    #"Stitch] Processing": lambda x: x.strip().split()[0].replace('[',''), 
    "[M::main] Real time": {"name": "minimap2", 'f': lambda x, struct: x.strip().split()[3]},

    "Predict]": {
        'f': lambda x, struct: catch_window(x.strip().split()[0].replace('[',''), struct),
        "f_post": lambda struct: post_window(struct),
    },
    "Sampler]": {
        'f': lambda x, struct: catch_window(x.strip().split()[0].replace('[',''), struct),
        "struct": "Predict]",
    },
    "Stitch]": {
        'f': lambda x, struct: catch_window(x.strip().split()[0].replace('[',''), struct),
        "f_post": lambda struct: post_window(struct),
    },
}

defaults = {
    "Predict]": {
        "first":None, "last":None, "days":0, "name":"predict", "wrapped":False
    },
    "Stitch]": {
        "first":None, "last":None, "days":0, "name":"stitch", "wrapped":False
    },
}


from datetime import datetime, timedelta
def post_window(struct):
    if struct["last"] < struct["first"]:
        struct["last"] = struct["last"] + timedelta(days=1)
    return (abs(struct["last"] - struct["first"]) + timedelta(days=struct["days"])).total_seconds()

def catch_window(x, struct):
    x = datetime.strptime(x, "%H:%M:%S")
    if not struct["first"]:
        struct["first"] = x

    if struct["last"]:
        if x < struct["last"]:
            struct["wrapped"] = True
        elif x > struct["first"] and struct["wrapped"]:
            # Next day
            struct["days"] += 1
            struct["wrapped"] = False
    struct["last"] = x
    return True


for log in os.scandir(LOG_DIR):

    if log.is_dir():
        continue
    res = {k:None for k in search}
    for query in search:
        if query in defaults:
            search[query].update(defaults[query])
    for line in open(os.path.join(LOG_DIR, log.name)):
        fields = line.strip().split()
        for query in search:
            if query in line:
                res[query] = search[query]["f"](line, search[ search[query].get("struct", query) ])

    broadstage = "default"
    if "pilon" in log.name:
        broadstage = "pilon"
    elif "medaka" in log.name:
        broadstage = "medaka"
    elif "racon" in log.name:
        broadstage = "racon"
    for query in res:
        if res[query]:
            if search[query].get("f_post", None):
                res[query] = search[query]["f_post"](search[query])
            print("\t".join([log.name, broadstage, search[query].get("name", query), str(res[query])]))

