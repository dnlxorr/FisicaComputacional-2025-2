N = 1.0
n_puntos = 100
dx = N / n_puntos
dy = N / n_puntos

x_val = [dx * i for i in range(n_puntos)]
y_val = [dy * i for i in range(n_puntos)]

# Mallas
X = [[x_val[i] for i in range(n_puntos)] for j in range(n_puntos)]
Y = [[y_val[j] for i in range(n_puntos)] for j in range(n_puntos)]

for y_val in Y:
    print(y_val)

#for y_val in Y:
#    print(y_val)