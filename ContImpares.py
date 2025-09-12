a = input('Ingrese un número cualquiera:')

def cont_impares(n):
    try:
        m = int(n)
        if m >= 0:
            c, i = 0, 1
            for i in range(m + 1):
                    if i % 2 == 1:
                        c += 1

            print(f'Hasta el número {n}, hay {c} números impares')
        else:
            print('Error: No se aceptan números negativos')
    except:
        return print('Error: No se aceptan letras, símbolos, decimales ni negativos')

cont_impares(a)