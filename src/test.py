import sys

nodes = [1,2,3,4]
values = map(lambda x : 4 ** x, range(len(nodes)))

def index2assignment(ind) :
    index = ind
    tmp = [0] * len(values)
    for i in sorted(range(len(values)), reverse=True) :
        tmp[i] = index / values[i]
        index = index % values[i]
    return tmp

def assignment2index(a) :
    tmp = 0
    for i in range(len(a)) :
        tmp += (a[i] * values[i])
    return tmp


x = eval(sys.argv[1])

index = assignment2index(x)
newx = index2assignment(index)

print "%s\n%s\n%s" % (x, index, newx)

print values
