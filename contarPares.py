numero = input("Ingrese un numero entero")

def cont_pares(nEntrada):
    try:

        numero = int(nEntrada)

        if numero <= 0:
            print("Error: ingrese un número entero positivo")
            return 0

        contador = 0
        for i in range(1, numero + 1):
            if i % 2 == 0:
                contador += 1

        print("Cantidad de números pares entre 1 y {numero}: {contador}")
        return contador

    except ValueError:
        print("Error: la entrada no es un número entero válido, no se aceptan letras, simbolos ni decimales. Ingresa"
              " un numero entero")
        return 0

cont_pares(numero)



