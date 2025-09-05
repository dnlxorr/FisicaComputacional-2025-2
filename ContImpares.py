def cont_impares(n):
    try:
        m = int(n)
        c, i = 0, 1
        for i in range(m):
            if i % 2 == 1:
                c += 1
        return print(c)
    except ValueError:
        return print('Error: no se ingresó un entero.')

cont_impares('hola')