from numpy import array, size

def derivate(x, y):
    result = []
    for i in range(size(x)-1):
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        result.append(dy/dx)

    return [x[1:], array(result)]

