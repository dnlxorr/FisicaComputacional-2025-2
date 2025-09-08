# contar números primos de manera sencilla

def primo(num):
    if num < 2:
        return False
    for i in range(2, num):
        if num % i == 0:
            return False
    else:
        return True

# pedir el número
n = input("Escriba un número entero: ")

# validar que no tenga letras ni decimales
if n.isdigit():
    n = int(n)
    print("Los números primos hasta", n, "son:")
    cont = 0
    for i in range(2, n+1):
        if primo(i):
            print(i)
            cont = cont + 1
    print("En total hay", cont, "números primos")
else:
    print("Solo puede escribir números enteros, nada de letras ni decimales.")
