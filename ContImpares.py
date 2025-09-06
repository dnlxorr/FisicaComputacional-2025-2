def cont_impares(n):
    try:
        if int(n):
            c, i = 0, 1
            for i in range(n):
                if i % 2 == 1:
                    c += 1
        print(c)
    except:
        return print('Error: No se aceptan letras, símbolos o decimales')

cont_impares(3.99)