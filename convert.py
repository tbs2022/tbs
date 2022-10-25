import pickle
f = open("tmppartition")

workers = f.readlines()

rt = {}

for i in range(len(workers)):
    rt[i] = workers[i]

pickle.dump(rt, open("map_for_maxinet", "wb"))
