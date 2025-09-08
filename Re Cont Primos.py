def es_primo(n):
    if n < 2: return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0: return False
    return True

while True:
    try:
        n = float(input("Límite para contar primos: "))
        if not n.is_integer(): print("Error: Solo enteros."); continue
        n = int(n)
        if n < 1: print("Error: Entero > 0."); continue
        break
    except ValueError: print("Error: Solo números.")

contador = sum(1 for x in range(1, n + 1) if es_primo(x))
print(f"Primos hasta {n}: {contador}")