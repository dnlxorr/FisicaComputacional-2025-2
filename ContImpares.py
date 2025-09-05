def cont_impares(n):
    c,i = 0,1
    for i in range(n):
        if i % 2 == 1:
            c += 1
    print(c)

try:
    cont_impares(2)
except:
    print('Error: No se aceptan letras o símbolos')